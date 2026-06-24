/*
    SPEAR - Sorting Petabytes of Environmental and Ancient Reads
    Copyright (C) 2026 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@sund.ku.dk>
    University of Copenhagen, Globe Institute, Section for GeoGenetics
    Oster Voldgade 5-7, 1350 Copenhagen K, Denmark
*/

#include "spear/align/bam_writer.hpp"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <string>

// htslib has its own extern "C" guards; no extra wrapping needed here.
#include <htslib/sam.h>

namespace spear::align {

// Verify our AlignmentFlags constants match the htslib BAM_F* values.
static_assert(
    AlignmentFlags::Paired        == BAM_FPAIRED,        "AlignmentFlags::Paired mismatch"
);
static_assert(
    AlignmentFlags::ProperPair    == BAM_FPROPER_PAIR,   "AlignmentFlags::ProperPair mismatch"
);
static_assert(
    AlignmentFlags::Unmapped      == BAM_FUNMAP,         "AlignmentFlags::Unmapped mismatch"
);
static_assert(
    AlignmentFlags::MateUnmapped  == BAM_FMUNMAP,        "AlignmentFlags::MateUnmapped mismatch"
);
static_assert(
    AlignmentFlags::Reverse       == BAM_FREVERSE,       "AlignmentFlags::Reverse mismatch"
);
static_assert(
    AlignmentFlags::MateReverse   == BAM_FMREVERSE,      "AlignmentFlags::MateReverse mismatch"
);
static_assert(
    AlignmentFlags::Read1         == BAM_FREAD1,         "AlignmentFlags::Read1 mismatch"
);
static_assert(
    AlignmentFlags::Read2         == BAM_FREAD2,         "AlignmentFlags::Read2 mismatch"
);
static_assert(
    AlignmentFlags::Secondary     == BAM_FSECONDARY,     "AlignmentFlags::Secondary mismatch"
);
static_assert(
    AlignmentFlags::QcFail        == BAM_FQCFAIL,        "AlignmentFlags::QcFail mismatch"
);
static_assert(
    AlignmentFlags::Duplicate     == BAM_FDUP,           "AlignmentFlags::Duplicate mismatch"
);
static_assert(
    AlignmentFlags::Supplementary == BAM_FSUPPLEMENTARY, "AlignmentFlags::Supplementary mismatch"
);

// =================================================================================================
//     BamHeader
// =================================================================================================

BamHeader::BamHeader( genesis::sequence::SequenceSet const& refs )
{
    references.reserve( refs.size() );
    for( auto const& seq : refs ) {
        references.emplace_back( seq.label(), static_cast<int64_t>( seq.length() ));
    }
}

// =================================================================================================
//     Helpers
// =================================================================================================

// Build the XA:Z: tag value from a list of secondary hits.
// Format per entry: rname,strand+pos,CIGAR,NM;
// pos is 1-based; sign encodes strand ('+' = forward, '-' = reverse).
static std::string build_xa_string(
    std::span<AlignmentHit const* const> hits,
    std::vector<std::pair<std::string, int64_t>> const& references
) {
    std::string xa;
    for( auto const* hit : hits ) {
        if( hit->ref_id >= 0 && hit->ref_id < static_cast<int32_t>( references.size() )) {
            xa += references[ static_cast<std::size_t>( hit->ref_id ) ].first;
        } else {
            xa += '*';
        }
        xa += ',';

        // 1-based position, signed by strand
        bool const reverse = ( hit->flags & AlignmentFlags::Reverse ) != 0;
        xa += ( reverse ? '-' : '+' );
        xa += std::to_string( hit->ref_pos + 1 );
        xa += ',';

        xa += cigar_to_string( hit->cigar );
        xa += ',';

        xa += std::to_string( hit->tag_nm.value_or( 0 ));
        xa += ';';
    }
    return xa;
}

// Apply HitSelectionPolicy and return ordered pointers into record.hits.
static std::vector<AlignmentHit const*> select_hits(
    ReadRecord const& record,
    HitSelectionPolicy const& policy
) {
    if( record.hits.empty() ) {
        return {};
    }

    // Find best AS for score_window filter.
    int best_as = record.hits[0].tag_as;
    if( policy.score_window >= 0 ) {
        for( auto const& h : record.hits ) {
            best_as = std::max( best_as, h.tag_as );
        }
    }

    std::vector<AlignmentHit const*> selected;
    selected.reserve( record.hits.size() );
    for( auto const& hit : record.hits ) {
        if( policy.min_mapq > 0 && hit.mapq < static_cast<uint8_t>( policy.min_mapq )) {
            continue;
        }
        if( policy.score_window >= 0 && hit.tag_as < best_as - policy.score_window ) {
            continue;
        }
        selected.push_back( &hit );
        if( policy.max_hits > 0 && static_cast<int>( selected.size() ) >= policy.max_hits ) {
            break;
        }
    }
    return selected;
}

// =================================================================================================
//     Impl
// =================================================================================================

struct BamWriter::Impl
{
    // Internal htslib objects for writing BAM/SAM/CRAM.
    samFile*          fp  = nullptr;
    sam_hdr_t*        hdr = nullptr;
    bam1_t*           rec = nullptr;

    // Original header and settings for the writer.
    BamHeader         header_;
    BamWriterSettings settings_;

    Impl( std::string const& path, BamHeader header, BamWriterSettings const& settings )
        : header_( std::move( header ) )
        , settings_( settings )
    {
        // Build the htslib mode string from format setting.
        std::string mode;
        switch( settings.format ) {
            case BamWriterSettings::Format::Sam:
                mode = "w";
                break;
            case BamWriterSettings::Format::Bam: {
                int const lvl = std::min( std::max( settings.bam_compression, 0 ), 9 );
                mode = "wb" + std::to_string( lvl );
                break;
            }
            case BamWriterSettings::Format::Cram:
                throw std::runtime_error( "BamWriter: CRAM output is not yet implemented" );
        }

        fp = sam_open( path.c_str(), mode.c_str() );
        if( !fp ) {
            throw std::runtime_error( "BamWriter: failed to open output file: " + path );
        }

        hdr = sam_hdr_init();
        if( !hdr ) {
            sam_close( fp );
            throw std::runtime_error( "BamWriter: failed to initialise SAM header" );
        }

        if( sam_hdr_add_line( hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL ) < 0 ) {
            throw std::runtime_error( "BamWriter: failed to add @HD line" );
        }
        for( auto const& [name, len] : header_.references ) {
            auto const len_str = std::to_string( len );
            if( sam_hdr_add_line( hdr, "SQ", "SN", name.c_str(), "LN", len_str.c_str(), NULL ) < 0 ) {
                throw std::runtime_error( "BamWriter: failed to add @SQ line for: " + name );
            }
        }
        if( header_.read_group ) {
            BamHeader::ReadGroup const& rg = *header_.read_group;
            int rg_ret;
            if( rg.library.empty() ) {
                rg_ret = sam_hdr_add_line( hdr, "RG",
                    "ID", rg.id.c_str(), "SM", rg.sample.c_str(),
                    "PL", rg.platform.c_str(), NULL
                );
            } else {
                rg_ret = sam_hdr_add_line( hdr, "RG",
                    "ID", rg.id.c_str(), "SM", rg.sample.c_str(),
                    "LB", rg.library.c_str(), "PL", rg.platform.c_str(), NULL
                );
            }
            if( rg_ret < 0 ) {
                throw std::runtime_error( "BamWriter: failed to add @RG line for ID: " + rg.id );
            }
        }
        if( sam_hdr_add_line( hdr, "PG", "ID", "spear", "PN", "spear", NULL ) < 0 ) {
            throw std::runtime_error( "BamWriter: failed to add @PG line" );
        }

        if( sam_hdr_write( fp, hdr ) < 0 ) {
            throw std::runtime_error( "BamWriter: failed to write SAM header to: " + path );
        }

        rec = bam_init1();
        if( !rec ) {
            throw std::runtime_error( "BamWriter: failed to allocate BAM record buffer" );
        }
    }

    ~Impl()
    {
        if( rec ) bam_destroy1( rec );
        if( hdr ) sam_hdr_destroy( hdr );
        if( fp  ) sam_close( fp );
    }

    Impl( Impl const& ) = delete;
    Impl& operator=( Impl const& ) = delete;

    // Fill `rec` with data from `rr` and `hit` using the given `flags`.
    // Does NOT write to the file; call flush_record() afterwards.
    void fill_record( ReadRecord const& rr, AlignmentHit const& hit, uint16_t flags )
    {
        auto const& name = rr.read.label();
        auto const& seq  = rr.read.sites();
        auto const& qual = rr.read.phred_scores();

        int const ret = bam_set1(
            rec,
            name.size() + 1,
            name.c_str(),
            flags,
            hit.ref_id,
            static_cast<hts_pos_t>( hit.ref_pos ),
            hit.mapq,
            hit.cigar.size(),
            hit.cigar.empty() ? nullptr : hit.cigar.data(),
            hit.mate_ref_id,
            static_cast<hts_pos_t>( hit.mate_ref_pos ),
            static_cast<hts_pos_t>( hit.template_len ),
            seq.size(),
            seq.c_str(),
            qual.empty() ? nullptr : reinterpret_cast<char const*>( qual.data() ),
            0
        );
        if( ret < 0 ) {
            throw std::runtime_error( "BamWriter: bam_set1 failed for read: " + name );
        }

        bam_aux_update_int( rec, "AS", static_cast<int64_t>( hit.tag_as ));
        if( hit.tag_nm ) {
            bam_aux_update_int( rec, "NM", static_cast<int64_t>( *hit.tag_nm ));
        }
        if( hit.tag_md ) {
            bam_aux_update_str( rec, "MD", -1, hit.tag_md->c_str() );
        }
        if( rr.tag_rg ) {
            bam_aux_update_str( rec, "RG", -1, rr.tag_rg->c_str() );
        }
    }

    void flush_record( std::string const& read_name )
    {
        if( sam_write1( fp, hdr, rec ) < 0 ) {
            throw std::runtime_error( "BamWriter: sam_write1 failed for read: " + read_name );
        }
    }

    void write( ReadRecord&& rr )
    {
        auto const selected = select_hits( rr, settings_.hit_selection );
        if( selected.empty() ) {
            return;
        }

        if( settings_.multi_hit_mode == BamWriterSettings::MultiHitMode::SeparateRecords ) {
            // First hit is primary; all subsequent hits get AlignmentFlags::Secondary added.
            fill_record( rr, *selected[0], selected[0]->flags );
            flush_record( rr.read.label() );
            for( std::size_t i = 1; i < selected.size(); ++i ) {
                fill_record( rr, *selected[i], selected[i]->flags | AlignmentFlags::Secondary );
                flush_record( rr.read.label() );
            }
        } else if( settings_.multi_hit_mode == BamWriterSettings::MultiHitMode::XaTag ) {
            // XaTag: fill primary, append XA tag with secondaries, then write once.
            fill_record( rr, *selected[0], selected[0]->flags );
            if( selected.size() > 1 ) {
                std::string const xa = build_xa_string(
                    std::span( selected ).subspan( 1 ), header_.references
                );
                bam_aux_update_str( rec, "XA", -1, xa.c_str() );
            }
            flush_record( rr.read.label() );
        } else {
            throw std::invalid_argument( "BamWriter: unhandled MultiHitMode value" );
        }
    }
};

// =================================================================================================
//     BamWriter
// =================================================================================================

BamWriter::BamWriter(
    std::string const& path,
    BamHeader header,
    BamWriterSettings settings
)
    : impl_( std::make_unique<Impl>( path, std::move( header ), settings ))
{}

BamWriter::~BamWriter() = default;

BamWriter::BamWriter( BamWriter&& ) noexcept = default;
BamWriter& BamWriter::operator=( BamWriter&& ) noexcept = default;

void BamWriter::write( ReadRecord&& record )
{
    impl_->write( std::move( record ));
}

} // namespace spear::align
