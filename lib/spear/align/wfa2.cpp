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

#include "spear/align/wfa2.hpp"

#include <stdexcept>

extern "C" {
#include "wfa2/wavefront/wfa.h"
#include "wfa2/alignment/cigar.h"
}

namespace spear::align {

// =================================================================================================
//     Damage match function
// =================================================================================================

// Args passed into the WFA2 lambda: sequence data and damage zone boundaries.
struct DamageArgs
{
    char const* query;   // pattern (read) sequence
    char const* target;  // text (reference) sequence
    int ct_end;          // positions [0, ct_end): C→T treated as free match (5' damage)
    int ga_start;        // positions [ga_start, query_len): G→A treated as free match (3' damage)
};

// Branchless damage-aware match function for wavefront_align_lambda().
// Returns 1 (match) for exact matches plus tolerated damage substitutions within their zones.
static int damage_match_fn( int v, int h, void* args_ )
{
    auto const* a = static_cast<DamageArgs const*>( args_ );
    auto const  q = static_cast<unsigned char>( a->query[v]  );
    auto const  r = static_cast<unsigned char>( a->target[h] );
    return (q == r)
         | (( q == 'T' ) & ( r == 'C' ) & ( v <  a->ct_end   ))
         | (( q == 'A' ) & ( r == 'G' ) & ( v >= a->ga_start ));
}

// =================================================================================================
//     Impl
// =================================================================================================

struct Wfa2Aligner::Impl
{
    Wfa2Settings         settings;
    bool                 damage_enabled;
    wavefront_aligner_t* wfa = nullptr;

    explicit Impl( Wfa2Settings s )
        : settings( s )
        , damage_enabled( s.damage_ct_end != 0 || s.damage_ga_start != 0 )
    {
        // Validity check
        if( s.mismatch < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::mismatch must be >= 0" );
        }
        if( s.gap_open < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::gap_open must be >= 0" );
        }
        if( s.gap_extend < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::gap_extend must be >= 0" );
        }
        if( s.damage_ct_end < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::damage_ct_end must be >= 0" );
        }
        if( s.damage_ga_start < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::damage_ga_start must be >= 0" );
        }
        if( s.max_steps < 0 ) {
            throw std::invalid_argument( "Wfa2Settings::max_steps must be >= 0" );
        }

        // Transfer settings to WFA2's C API
        wavefront_aligner_attr_t attrs = wavefront_aligner_attr_default;
        attrs.distance_metric  = gap_affine;
        attrs.alignment_scope  = ( s.scope == Wfa2Settings::Scope::Score )
                                 ? compute_score : compute_alignment;
        attrs.alignment_form.span               = alignment_endsfree;
        attrs.alignment_form.pattern_begin_free = 0;
        attrs.alignment_form.pattern_end_free   = 0;
        attrs.alignment_form.text_begin_free    = 0;
        attrs.alignment_form.text_end_free      = 0;
        attrs.affine_penalties.match            = 0;
        attrs.affine_penalties.mismatch         = s.mismatch;
        attrs.affine_penalties.gap_opening      = s.gap_open;
        attrs.affine_penalties.gap_extension    = s.gap_extend;
        attrs.system.verbose                    = 0;
        attrs.system.max_alignment_steps        = s.max_steps;

        // Create the WFA2 aligner
        wfa = wavefront_aligner_new( &attrs );
        wavefront_aligner_set_heuristic_none( wfa );
    }

    ~Impl()
    {
        if( wfa ) {
            wavefront_aligner_delete( wfa );
        }
    }

    Impl( Impl const& ) = delete;
    Impl& operator=( Impl const& ) = delete;

    inline int run_align( std::string const& query, std::string const& target )
    {
        // Set free-end alignment bounds and call the appropriate WFA2 functions.
        // With damage at the end, we use the lambda variant to pass in the match function.
        // Otherwise, we call the slightly faster standard align function.
        int const qlen = static_cast<int>( query.size()  );
        int const tlen = static_cast<int>( target.size() );
        wavefront_aligner_set_alignment_free_ends( wfa, 0, 0, tlen, tlen );
        if( damage_enabled ) {
            DamageArgs args{
                query.c_str(), target.c_str(),
                settings.damage_ct_end,
                settings.damage_ga_start
            };
            return wavefront_align_lambda( wfa, damage_match_fn, &args, qlen, tlen );
        }
        return wavefront_align( wfa, query.c_str(), qlen, target.c_str(), tlen );
    }
};

// =================================================================================================
//     Helpers
// =================================================================================================

static Wfa2Result::Status map_wfa_status( int wfa_status )
{
    switch( wfa_status ) {
        case WF_STATUS_ALG_COMPLETED:     return Wfa2Result::Status::Ok;
        case WF_STATUS_MAX_STEPS_REACHED: return Wfa2Result::Status::ScoreCapExceeded;
        case WF_STATUS_OOM:               return Wfa2Result::Status::OutOfMemory;
        default:                          return Wfa2Result::Status::Failed;
    }
}

// =================================================================================================
//     Wfa2Aligner
// =================================================================================================

Wfa2Aligner::Wfa2Aligner( Wfa2Settings settings )
    : impl_( std::make_unique<Impl>( settings ))
{}

Wfa2Aligner::~Wfa2Aligner() = default;

Wfa2Aligner::Wfa2Aligner( Wfa2Aligner&& ) noexcept = default;
Wfa2Aligner& Wfa2Aligner::operator=( Wfa2Aligner&& ) noexcept = default;

// -------------------------------------------------------------------------

Wfa2Result Wfa2Aligner::align_score( std::string const& query, std::string const& target )
{
    // Sanity check
    if( impl_->settings.scope != Wfa2Settings::Scope::Score ) {
        throw std::logic_error( "align_score() called on a Scope::Cigar aligner" );
    }

    // Compute alignment
    auto const s = map_wfa_status( impl_->run_align( query, target ));
    if( s != Wfa2Result::Status::Ok ) {
        return Wfa2Result{ .status = s };
    }

    // Copy results to our struct
    auto* wfa = impl_->wfa;
    return Wfa2Result{
        .ref_begin = 0,                              // not computed in score mode
        .ref_end   = wfa->alignment_end_pos.offset,  // fitting endpoint in reference
        .score     = wfa->cigar->score,
        .status    = Wfa2Result::Status::Ok,
    };
}

Wfa2Result Wfa2Aligner::align_cigar( std::string const& query, std::string const& target )
{
    // Sanity check
    if( impl_->settings.scope != Wfa2Settings::Scope::Cigar ) {
        throw std::logic_error( "align_cigar() called on a Scope::Score aligner" );
    }

    // Compute alignment
    auto const s = map_wfa_status( impl_->run_align( query, target ));
    if( s != Wfa2Result::Status::Ok ) {
        return Wfa2Result{ .status = s };
    }

    auto* wfa = impl_->wfa;

    // Encode into WFA2's internal buffer (BAM uint32_t RLE format) and get a pointer to it.
    uint32_t* cigar_buffer = nullptr;
    int       cigar_length = 0;
    cigar_get_CIGAR( wfa->cigar, false, &cigar_buffer, &cigar_length );

    // Derive ref_begin by counting reference-consuming ops (M, X, D).
    cigar_t const* cig = wfa->cigar;
    int ref_span = 0;
    for( int i = cig->begin_offset; i < cig->end_offset; ++i ) {
        char const op = cig->operations[i];
        if( op == 'M' || op == 'X' || op == 'D' ) {
            ++ref_span;
        }
    }

    // Copy results to our struct
    return Wfa2Result{
        .ref_begin = cig->end_h - ref_span,
        .ref_end   = cig->end_h,
        .score     = cig->score,
        .status    = Wfa2Result::Status::Ok,
        .cigar     = { cigar_buffer, static_cast<std::size_t>( cigar_length ) },
    };
}

} // namespace spear::align
