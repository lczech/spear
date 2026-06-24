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

/**
 * @brief Tests for BamWriter, BamHeader, and HitSelectionPolicy (bam_writer.hpp / bam_writer.cpp).
 *
 * All file-I/O tests write SAM format (Format::Sam) so records can be verified by reading the
 * plain-text output directly, without a separate htslib read-back pass.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/align/bam_writer.hpp"
#include "spear/align/read_record.hpp"

#include <genesis/sequence/sequence.hpp>
#include <genesis/sequence/sequence_set.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace spear::align;
using genesis::sequence::Sequence;
using genesis::sequence::SequenceSet;

// BAM CIGAR encoding helper: (length << 4) | CigarOp::Op
static constexpr uint32_t op( uint32_t len, uint32_t code ) { return (len << 4) | code; }

// =================================================================================================
//     Test helpers
// =================================================================================================

// Build a BamHeader from a flat list of (name, length) pairs.
static BamHeader make_header( std::vector<std::pair<std::string, int64_t>> refs )
{
    BamHeader hdr;
    hdr.references = std::move( refs );
    return hdr;
}

// Build a ReadRecord with a single AlignmentHit.
static ReadRecord make_record(
    std::string const&    name,
    std::string const&    seq,
    int32_t               ref_id  = 0,
    int64_t               ref_pos = 0,
    uint8_t               mapq    = 60,
    std::vector<uint32_t> cigar   = { op(4, CigarOp::E) },
    uint16_t              flags   = 0,
    int                   tag_as  = 0
) {
    ReadRecord rr;
    rr.read.label( name );
    rr.read.sites( seq );
    AlignmentHit hit;
    hit.ref_id  = ref_id;
    hit.ref_pos = ref_pos;
    hit.mapq    = mapq;
    hit.cigar   = std::move( cigar );
    hit.flags   = flags;
    hit.tag_as  = tag_as;
    rr.hits.push_back( std::move( hit ) );
    return rr;
}

// Return BamWriterSettings configured for SAM (plain-text) output.
static BamWriterSettings sam_settings()
{
    BamWriterSettings s;
    s.format = BamWriterSettings::Format::Sam;
    return s;
}

// Return all @-prefixed header lines from a SAM file.
static std::vector<std::string> read_header_lines( std::string const& path )
{
    std::ifstream f( path );
    std::vector<std::string> lines;
    std::string line;
    while( std::getline( f, line ) ) {
        if( !line.empty() && line[0] == '@' ) {
            lines.push_back( line );
        }
    }
    return lines;
}

// Return alignment records from a SAM file as vectors of tab-split fields.
static std::vector<std::vector<std::string>> read_sam_records( std::string const& path )
{
    std::ifstream f( path );
    std::vector<std::vector<std::string>> recs;
    std::string line;
    while( std::getline( f, line ) ) {
        if( line.empty() || line[0] == '@' ) continue;
        std::vector<std::string> fields;
        std::istringstream iss( line );
        std::string field;
        while( std::getline( iss, field, '\t' ) ) {
            fields.push_back( std::move( field ) );
        }
        recs.push_back( std::move( fields ) );
    }
    return recs;
}

// Find an aux field in fields[11:] by its 2-char tag name; returns "" if absent.
static std::string find_aux( std::vector<std::string> const& fields, std::string const& tag )
{
    for( std::size_t i = 11; i < fields.size(); ++i ) {
        if( fields[i].size() >= 2 && fields[i].substr( 0, 2 ) == tag ) {
            return fields[i];
        }
    }
    return {};
}

// Extract the value portion (after the last ':') from a TAG:TYPE:VALUE aux string.
static std::string aux_value( std::string const& aux_field )
{
    auto const pos = aux_field.rfind( ':' );
    if( pos == std::string::npos ) return {};
    return aux_field.substr( pos + 1 );
}

// Return true if any string in lines contains needle as a substring.
static bool header_has( std::vector<std::string> const& lines, std::string const& needle )
{
    for( auto const& l : lines ) {
        if( l.find( needle ) != std::string::npos ) return true;
    }
    return false;
}

// =================================================================================================
//     BamHeader
// =================================================================================================

TEST( BamHeaderTest, DefaultConstruct )
{
    BamHeader hdr;
    EXPECT_TRUE( hdr.references.empty() );
    EXPECT_FALSE( hdr.read_group.has_value() );
}

TEST( BamHeaderTest, FromSequenceSet )
{
    SequenceSet refs;
    refs.add( Sequence( "chr1", "ACGT" ) );
    refs.add( Sequence( "chr2", "ACGTACGT" ) );
    BamHeader hdr( refs );
    ASSERT_EQ( hdr.references.size(), 2 );
    EXPECT_EQ( hdr.references[0].first,  "chr1" );
    EXPECT_EQ( hdr.references[0].second, 4 );
    EXPECT_EQ( hdr.references[1].first,  "chr2" );
    EXPECT_EQ( hdr.references[1].second, 8 );
}

TEST( BamHeaderTest, EmptySequenceSet )
{
    SequenceSet refs;
    BamHeader hdr( refs );
    EXPECT_TRUE( hdr.references.empty() );
}

TEST( BamHeaderTest, ReadGroupFields )
{
    BamHeader hdr;
    hdr.read_group = BamHeader::ReadGroup{ "RG1", "Sample1", "Lib1", "ILLUMINA" };
    ASSERT_TRUE( hdr.read_group.has_value() );
    EXPECT_EQ( hdr.read_group->id,       "RG1" );
    EXPECT_EQ( hdr.read_group->sample,   "Sample1" );
    EXPECT_EQ( hdr.read_group->library,  "Lib1" );
    EXPECT_EQ( hdr.read_group->platform, "ILLUMINA" );
}

TEST( BamHeaderTest, ReadGroupDefaultPlatform )
{
    BamHeader::ReadGroup rg;
    EXPECT_EQ( rg.platform, "ILLUMINA" );
}

// =================================================================================================
//     BamWriter construction
// =================================================================================================

TEST( BamWriterConstruct, OpenSam )
{
    TempFile tmp;
    BamHeader hdr = make_header({{"chr1", 100}});
    EXPECT_NO_THROW( BamWriter w( tmp.path, std::move(hdr), sam_settings() ));
}

TEST( BamWriterConstruct, OpenBam )
{
    TempFile tmp;
    BamHeader hdr = make_header({{"chr1", 100}});
    EXPECT_NO_THROW( BamWriter w( tmp.path, std::move(hdr) ));
}

TEST( BamWriterConstruct, CramThrows )
{
    TempFile tmp;
    BamHeader hdr = make_header({{"chr1", 100}});
    BamWriterSettings s;
    s.format = BamWriterSettings::Format::Cram;
    EXPECT_THROW( BamWriter( tmp.path, std::move(hdr), s ), std::runtime_error );
}

TEST( BamWriterConstruct, BadPathThrows )
{
    BamHeader hdr = make_header({{"chr1", 100}});
    EXPECT_THROW(
        BamWriter( "/nonexistent_spear_test_dir/out.sam", std::move(hdr), sam_settings() ),
        std::runtime_error
    );
}

TEST( BamWriterConstruct, MoveConstruct )
{
    TempFile tmp;
    BamHeader hdr = make_header({{"chr1", 100}});
    BamWriter w1( tmp.path, std::move(hdr), sam_settings() );
    BamWriter w2( std::move(w1) );
    // w2 should be fully usable after the move
    EXPECT_NO_THROW( w2.write( make_record("r1", "ACGT") ));
}

// =================================================================================================
//     SAM header content
// =================================================================================================

TEST( BamWriterHeader, HdSqPgLines )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}, {"chr2", 200}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
    }
    auto lines = read_header_lines( tmp.path );
    EXPECT_TRUE( header_has( lines, "@HD" ));
    EXPECT_TRUE( header_has( lines, "VN:1.6" ));
    EXPECT_TRUE( header_has( lines, "SO:unsorted" ));
    EXPECT_TRUE( header_has( lines, "SN:chr1" ));
    EXPECT_TRUE( header_has( lines, "LN:100" ));
    EXPECT_TRUE( header_has( lines, "SN:chr2" ));
    EXPECT_TRUE( header_has( lines, "LN:200" ));
    EXPECT_TRUE( header_has( lines, "@PG" ));
    EXPECT_TRUE( header_has( lines, "ID:spear" ));
}

TEST( BamWriterHeader, RgLinePresent )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        hdr.read_group = BamHeader::ReadGroup{ "RG1", "Sample1", "", "ILLUMINA" };
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
    }
    auto lines = read_header_lines( tmp.path );
    EXPECT_TRUE( header_has( lines, "@RG" ));
    EXPECT_TRUE( header_has( lines, "ID:RG1" ));
    EXPECT_TRUE( header_has( lines, "SM:Sample1" ));
    EXPECT_TRUE( header_has( lines, "PL:ILLUMINA" ));
}

TEST( BamWriterHeader, RgLineWithLibrary )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        hdr.read_group = BamHeader::ReadGroup{ "RG1", "Sample1", "Lib1", "ILLUMINA" };
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
    }
    auto lines = read_header_lines( tmp.path );
    EXPECT_TRUE( header_has( lines, "LB:Lib1" ));
}

TEST( BamWriterHeader, RgLineAbsent )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
    }
    auto lines = read_header_lines( tmp.path );
    EXPECT_FALSE( header_has( lines, "@RG" ));
}

// =================================================================================================
//     HitSelectionPolicy
// =================================================================================================

TEST( HitSelection, EmptyHits )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        // no hits added
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 0 );
}

TEST( HitSelection, NoFiltersAllPass )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 3; ++i ) {
            AlignmentHit hit;
            hit.ref_id  = 0;
            hit.ref_pos = i * 10;
            hit.cigar   = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 3 );
}

TEST( HitSelection, MaxHitsLimit )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.max_hits = 2;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 4; ++i ) {
            AlignmentHit hit;
            hit.ref_id  = 0;
            hit.ref_pos = i * 10;
            hit.cigar   = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 2 );
}

TEST( HitSelection, MinMapqFiltersLowQuality )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.min_mapq = 30;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit low, high;
        low.ref_id  = 0; low.ref_pos  = 0;  low.mapq  = 25; low.cigar  = { op(4, CigarOp::E) };
        high.ref_id = 0; high.ref_pos = 10; high.mapq = 30; high.cigar = { op(4, CigarOp::E) };
        rr.hits.push_back( low );
        rr.hits.push_back( high );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][4], "30" ); // only the higher-mapq hit survived
}

TEST( HitSelection, MinMapqAtThresholdPasses )
{
    // The filter is `mapq < min_mapq`, so a hit with mapq == min_mapq must pass.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.min_mapq = 30;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit hit;
        hit.ref_id = 0; hit.ref_pos = 0; hit.mapq = 30; hit.cigar = { op(4, CigarOp::E) };
        rr.hits.push_back( hit );
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 1 );
}

TEST( HitSelection, MinMapqZeroDisablesFilter )
{
    // When min_mapq == 0 the filter is bypassed entirely (gate: `policy.min_mapq > 0`).
    // A hit with mapq=0 must still be written.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.min_mapq = 0;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit hit;
        hit.ref_id = 0; hit.ref_pos = 0; hit.mapq = 0; hit.cigar = { op(4, CigarOp::E) };
        rr.hits.push_back( hit );
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 1 );
}

TEST( HitSelection, ScoreWindowDropsDistantHits )
{
    // score_window=5: drop hits with AS < best_AS - 5.
    // Hits: AS=0 (best), AS=-4 (passes: -4 >= -5), AS=-6 (dropped: -6 < -5).
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.score_window = 5;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int as : { 0, -4, -6 } ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = 0; hit.cigar = { op(4, CigarOp::E) }; hit.tag_as = as;
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 2 );
}

TEST( HitSelection, ScoreWindowBoundaryPasses )
{
    // Filter is strictly `tag_as < best_as - score_window`, so the boundary value must pass.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.score_window = 5;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int as : { 0, -5 } ) { // -5 is exactly best_as - window: 0 - 5 = -5
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = 0; hit.cigar = { op(4, CigarOp::E) }; hit.tag_as = as;
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 2 );
}

TEST( HitSelection, ScoreWindowZeroKeepsOnlyBest )
{
    // score_window=0: only hits with AS == best_AS survive.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.score_window = 0;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int as : { 0, -1 } ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = 0; hit.cigar = { op(4, CigarOp::E) }; hit.tag_as = as;
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 1 );
}

TEST( HitSelection, ScoreWindowNegativeDisablesFilter )
{
    // score_window=-1 (default) disables the AS filter entirely.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.score_window = -1;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int as : { 0, -100, -9999 } ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = 0; hit.cigar = { op(4, CigarOp::E) }; hit.tag_as = as;
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 3 );
}

TEST( HitSelection, CombinedMaxHitsAndMinMapq )
{
    // max_hits=2, min_mapq=30: from 4 hits (mapq 20,35,40,50) only 3 pass mapq,
    // then max_hits caps at 2.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriterSettings s = sam_settings();
        s.hit_selection.max_hits = 2;
        s.hit_selection.min_mapq = 30;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( uint8_t mq : { 20, 35, 40, 50 } ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = 0; hit.mapq = mq; hit.cigar = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 2 );
}

// =================================================================================================
//     MultiHitMode::SeparateRecords
// =================================================================================================

TEST( SeparateRecords, SingleHitNoSecondaryFlag )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::SeparateRecords;
        BamWriter w( tmp.path, std::move(hdr), s );
        w.write( make_record("r1", "ACGT") );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( std::stoi( recs[0][1] ) & AlignmentFlags::Secondary, 0 );
}

TEST( SeparateRecords, TwoHitsSecondaryFlag )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::SeparateRecords;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 2; ++i ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = i * 10; hit.cigar = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 2 );
    EXPECT_EQ( std::stoi( recs[0][1] ) & AlignmentFlags::Secondary, 0 );
    EXPECT_NE( std::stoi( recs[1][1] ) & AlignmentFlags::Secondary, 0 );
}

TEST( SeparateRecords, ThreeHitsAllSecondariesFlagged )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::SeparateRecords;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 3; ++i ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = i * 10; hit.cigar = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 3 );
    EXPECT_EQ( std::stoi( recs[0][1] ) & AlignmentFlags::Secondary, 0 );
    EXPECT_NE( std::stoi( recs[1][1] ) & AlignmentFlags::Secondary, 0 );
    EXPECT_NE( std::stoi( recs[2][1] ) & AlignmentFlags::Secondary, 0 );
}

TEST( SeparateRecords, SecondaryPreservesOtherFlags )
{
    // Secondary hit already has Reverse; BamWriter should OR in Secondary without clearing Reverse.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::SeparateRecords;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit primary;
        primary.ref_id = 0; primary.ref_pos = 0; primary.cigar = { op(4, CigarOp::E) };
        AlignmentHit secondary;
        secondary.ref_id = 0; secondary.ref_pos = 10; secondary.cigar = { op(4, CigarOp::E) };
        secondary.flags = AlignmentFlags::Reverse;
        rr.hits.push_back( primary );
        rr.hits.push_back( secondary );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 2 );
    int flag1 = std::stoi( recs[1][1] );
    EXPECT_NE( flag1 & AlignmentFlags::Reverse,   0 );
    EXPECT_NE( flag1 & AlignmentFlags::Secondary, 0 );
}

// =================================================================================================
//     MultiHitMode::XaTag
// =================================================================================================

TEST( XaTag, SingleHitNoXaTag )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        w.write( make_record("r1", "ACGT") );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( find_aux( recs[0], "XA" ), "" );
}

TEST( XaTag, MultipleHitsOnlyOnePrimaryRecord )
{
    // XaTag mode emits exactly one BAM record regardless of hit count.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 3; ++i ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = i * 10; hit.cigar = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 1 );
}

TEST( XaTag, TwoHitsXaEntryFormat )
{
    // Secondary at ref_pos=20 (1-based=21), forward strand, NM=1 → "chr1,+21,4=,1;"
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit h1, h2;
        h1.ref_id = 0; h1.ref_pos = 5;  h1.cigar = { op(4, CigarOp::E) };
        h2.ref_id = 0; h2.ref_pos = 20; h2.cigar = { op(4, CigarOp::E) }; h2.tag_nm = 1;
        rr.hits.push_back( h1 );
        rr.hits.push_back( h2 );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    auto xa = find_aux( recs[0], "XA" );
    ASSERT_FALSE( xa.empty() );
    EXPECT_NE( aux_value( xa ).find("chr1,+21,4=,1;"), std::string::npos );
}

TEST( XaTag, MultipleSecondariesTwoEntries )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 1000}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        for( int i = 0; i < 3; ++i ) {
            AlignmentHit hit;
            hit.ref_id = 0; hit.ref_pos = i * 100; hit.cigar = { op(4, CigarOp::E) };
            rr.hits.push_back( hit );
        }
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    // Two secondary entries means the XA value contains exactly two ';' delimiters.
    auto xa_val = aux_value( find_aux( recs[0], "XA" ) );
    EXPECT_EQ( std::count( xa_val.begin(), xa_val.end(), ';' ), 2 );
}

TEST( XaTag, ReverseStrandPrefixedWithMinus )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit h1, h2;
        h1.ref_id = 0; h1.ref_pos = 0;  h1.cigar = { op(4, CigarOp::E) };
        h2.ref_id = 0; h2.ref_pos = 10; h2.cigar = { op(4, CigarOp::E) };
        h2.flags = AlignmentFlags::Reverse;
        rr.hits.push_back( h1 );
        rr.hits.push_back( h2 );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    // ref_pos=10 → 1-based=11; reverse strand → '-'
    EXPECT_NE( aux_value( find_aux( recs[0], "XA" )).find("chr1,-11,"), std::string::npos );
}

TEST( XaTag, NmFallbackToZeroWhenAbsent )
{
    // Secondaries without tag_nm should emit ,0; in XA, not an empty field.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit h1, h2;
        h1.ref_id = 0; h1.ref_pos = 0;  h1.cigar = { op(4, CigarOp::E) };
        h2.ref_id = 0; h2.ref_pos = 10; h2.cigar = { op(4, CigarOp::E) };
        // h2.tag_nm intentionally left as nullopt
        rr.hits.push_back( h1 );
        rr.hits.push_back( h2 );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_NE( aux_value( find_aux( recs[0], "XA" )).find(",0;"), std::string::npos );
}

TEST( XaTag, UnmappedSecondaryUsesStarRname )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriterSettings s = sam_settings();
        s.multi_hit_mode = BamWriterSettings::MultiHitMode::XaTag;
        BamWriter w( tmp.path, std::move(hdr), s );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit h1, h2;
        h1.ref_id = 0; h1.ref_pos = 0; h1.cigar = { op(4, CigarOp::E) };
        h2.ref_id = -1; // unmapped
        rr.hits.push_back( h1 );
        rr.hits.push_back( h2 );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_NE( aux_value( find_aux( recs[0], "XA" )).find("*,"), std::string::npos );
}

// =================================================================================================
//     Record fields
// =================================================================================================

TEST( RecordFields, QnameSeqCigar )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("myread", "ACGT", 0, 5, 60, { op(2, CigarOp::E), op(2, CigarOp::X) }) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][0], "myread" ); // QNAME
    EXPECT_EQ( recs[0][9], "ACGT" );   // SEQ
    EXPECT_EQ( recs[0][5], "2=2X" );   // CIGAR
}

TEST( RecordFields, RefnameAndPos )
{
    // ref_id=1 → "chr2"; ref_pos=49 → SAM POS 50 (1-based).
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}, {"chr2", 200}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT", 1, 49, 60, { op(4, CigarOp::E) }) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][2], "chr2" );
    EXPECT_EQ( recs[0][3], "50" );
}

TEST( RecordFields, FlagAndMapq )
{
    // Reverse flag (0x10 = 16) must appear verbatim in FLAG; mapq must be preserved.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT", 0, 0, 37, { op(4, CigarOp::E) }, AlignmentFlags::Reverse) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][1], "16" );
    EXPECT_EQ( recs[0][4], "37" );
}

TEST( RecordFields, QualPresentPhredPlusThirtyThree )
{
    // htslib converts raw Phred scores to Phred+33 ASCII when writing SAM.
    // Storing {40,40,40,40} → each byte 40+33=73='I' in the QUAL column.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr = make_record("r1", "ACGT");
        rr.read.phred_scores({ 40, 40, 40, 40 });
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][10], "IIII" );
}

TEST( RecordFields, QualAbsentWritesStar )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT") ); // no phred_scores
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][10], "*" );
}

TEST( RecordFields, AsTagAlwaysPresentWithCorrectValue )
{
    // AS is always written; value may be negative (WFA2 penalty convention: 0=perfect, <0=worse).
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT", 0, 0, 60, { op(4, CigarOp::E) }, 0, -7) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_NE( find_aux( recs[0], "AS" ), "" );
    EXPECT_EQ( aux_value( find_aux( recs[0], "AS" )), "-7" );
}

TEST( RecordFields, NmTagPresentWhenSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr = make_record("r1", "ACGT");
        rr.hits[0].tag_nm = 2;
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( aux_value( find_aux( recs[0], "NM" )), "2" );
}

TEST( RecordFields, NmTagAbsentWhenNotSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        // tag_nm is nullopt by default
        w.write( make_record("r1", "ACGT") );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( find_aux( recs[0], "NM" ), "" );
}

TEST( RecordFields, MdTagPresentWhenSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr = make_record("r1", "ACGT");
        rr.hits[0].tag_md = "3T0";
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( aux_value( find_aux( recs[0], "MD" )), "3T0" );
}

TEST( RecordFields, MdTagAbsentWhenNotSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT") );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( find_aux( recs[0], "MD" ), "" );
}

TEST( RecordFields, RgTagPresentWhenSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        hdr.read_group = BamHeader::ReadGroup{ "RG1", "S1", "", "ILLUMINA" };
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr = make_record("r1", "ACGT");
        rr.tag_rg = "RG1";
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( aux_value( find_aux( recs[0], "RG" )), "RG1" );
}

TEST( RecordFields, RgTagAbsentWhenNotSet )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT") );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( find_aux( recs[0], "RG" ), "" );
}

// =================================================================================================
//     Edge cases
// =================================================================================================

TEST( EdgeCase, UnmappedRead )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        ReadRecord rr;
        rr.read.label("r1");
        rr.read.sites("ACGT");
        AlignmentHit hit;
        hit.ref_id  = -1;
        hit.ref_pos = -1; // SAM POS=0 requires pos<0 in htslib (0-based -1 → output 0)
        hit.flags   = AlignmentFlags::Unmapped;
        // cigar is empty (default)
        rr.hits.push_back( hit );
        w.write( std::move(rr) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][2], "*" ); // RNAME
    EXPECT_EQ( recs[0][3], "0" ); // POS=0 for unmapped
    EXPECT_EQ( recs[0][5], "*" ); // CIGAR
}

TEST( EdgeCase, MultipleReferencesCorrectRname )
{
    // ref_id=2 → should map to the third reference name.
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 100}, {"chr2", 200}, {"chr3", 300}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        w.write( make_record("r1", "ACGT", 2, 0, 60, { op(4, CigarOp::E) }) );
    }
    auto recs = read_sam_records( tmp.path );
    ASSERT_EQ( recs.size(), 1 );
    EXPECT_EQ( recs[0][2], "chr3" );
}

TEST( EdgeCase, MultipleConsecutiveWrites )
{
    TempFile tmp;
    {
        BamHeader hdr = make_header({{"chr1", 10000}});
        BamWriter w( tmp.path, std::move(hdr), sam_settings() );
        for( int i = 0; i < 50; ++i ) {
            w.write( make_record("r" + std::to_string(i), "ACGT", 0, i * 10, 60, { op(4, CigarOp::E) }) );
        }
    }
    EXPECT_EQ( read_sam_records( tmp.path ).size(), 50 );
}
