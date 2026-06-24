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
#include "spear/align/cigar.hpp"

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
static inline int damage_match_fn( int v, int h, void* args_ )
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
    if( !query.empty() && target.empty() ) {
        throw std::invalid_argument( "wfa2::align_score: non-empty query against empty reference" );
    }
    if( query.empty() ) {
        return Wfa2Result{
            .ref_begin = 0, .ref_end = 0, .score = 0, .status = Wfa2Result::Status::Ok, .cigar = {}
        };
    }

    // Compute alignment
    auto const s = map_wfa_status( impl_->run_align( query, target ));
    if( s != Wfa2Result::Status::Ok ) {
        return Wfa2Result{ .status = s, .cigar = {} };
    }

    // Copy results to our struct
    auto* wfa = impl_->wfa;
    return Wfa2Result{
        .ref_begin = 0,                              // not computed in score mode
        .ref_end   = wfa->alignment_end_pos.offset,  // fitting endpoint in reference
        .score     = wfa->cigar->score,
        .status    = Wfa2Result::Status::Ok,
        .cigar     = {},
    };
}

Wfa2Result Wfa2Aligner::align_cigar( std::string const& query, std::string const& target )
{
    // Sanity check
    if( impl_->settings.scope != Wfa2Settings::Scope::Cigar ) {
        throw std::logic_error( "align_cigar() called on a Scope::Score aligner" );
    }
    if( !query.empty() && target.empty() ) {
        throw std::invalid_argument( "wfa2::align_cigar: non-empty query against empty reference" );
    }
    if( query.empty() ) {
        return Wfa2Result{
            .ref_begin = 0, .ref_end = 0, .score = 0, .status = Wfa2Result::Status::Ok, .cigar = {}
        };
    }

    // Compute alignment
    auto const s = map_wfa_status( impl_->run_align( query, target ));
    if( s != Wfa2Result::Status::Ok ) {
        return Wfa2Result{ .status = s, .cigar = {} };
    }

    auto* wfa = impl_->wfa;
    cigar_t* cig = wfa->cigar;

    // Encode into WFA2's internal buffer (BAM uint32_t RLE format) and get a pointer to it.
    // cigar_buffer points into cig->cigar_buffer (no allocation); copy into the vector below.
    uint32_t* cigar_buffer = nullptr;
    int       cigar_length = 0;
    cigar_get_CIGAR( cig, impl_->settings.use_extended_cigar, &cigar_buffer, &cigar_length );

    // WFA2 names its ops from the pattern-matching convention (transform pattern into text),
    // which is the opposite of SAM (transform reference→query). In WFA2's output:
    //   op 1 ("I") = text/reference consumes without pattern/query  (= SAM D)
    //   op 2 ("D") = pattern/query consumes without text/reference  (= SAM I)
    // In the glocal/fitting mode (text_begin_free=tlen, text_end_free=tlen), unaligned
    // reference positions at both flanks appear as leading/trailing op-1 runs.
    //
    // We correct this in one forward pass + back-trim:
    //   1. Skip leading op-1 runs → their total length is ref_begin.
    //   2. For the remaining ops, swap op 1↔2 to produce SAM-compliant I/D.
    //   3. Trim trailing op-2 (post-swap) runs from the output → they were the trailing clips.
    //   4. Sanity-check query_span == qlen and ref_end == cig->end_h - trailing_clip_len.
    //
    // WFA2 gap-affine only emits M(0), I(1), D(2), E(7), X(8); others indicate a bug.
    // Note: WFA2's I(1) and D(2) are swapped vs SAM convention; see Step 2 below.

    int const qlen = static_cast<int>( query.size() );
    int ref_begin  = 0;
    int ref_span   = 0;
    int query_span = 0;
    int cur_op_idx = 0;

    // Step 1: count leading WFA2-I ops (= SAM-D, free-end reference clips before alignment)
    while( cur_op_idx < cigar_length && (cigar_buffer[cur_op_idx] & 0xF) == CigarOp::I ) {
        ref_begin += static_cast<int>( cigar_buffer[cur_op_idx] >> 4 );
        ++cur_op_idx;
    }

    // Step 2: process remaining ops — swap I↔D, track ref/query span, copy to output vector
    std::vector<uint32_t> cigar;
    cigar.reserve( static_cast<size_t>( cigar_length - cur_op_idx ));
    for( ; cur_op_idx < cigar_length; ++cur_op_idx ) {
        uint32_t op  = cigar_buffer[cur_op_idx] & 0xF;
        uint32_t len = cigar_buffer[cur_op_idx] >> 4;

        // Swap WFA2 convention → SAM convention
        if( op == CigarOp::I ) {
            op = CigarOp::D;  // WFA2-I (ref-consuming) → SAM D
        } else if( op == CigarOp::D ) {
            op = CigarOp::I;  // WFA2-D (query-consuming) → SAM I
        } else if( op != CigarOp::M && op != CigarOp::E && op != CigarOp::X ) {
            throw std::logic_error( "wfa2: unexpected CIGAR op code " + std::to_string( op ));
        }

        // SAM: M, D, E, X consume reference; M, I, E, X consume query
        if( op == CigarOp::M || op == CigarOp::D || op == CigarOp::E || op == CigarOp::X ) {
            ref_span   += static_cast<int>( len );
        }
        if( op == CigarOp::M || op == CigarOp::I || op == CigarOp::E || op == CigarOp::X ) {
            query_span += static_cast<int>( len );
        }

        cigar.push_back( (len << 4) | op );
    }

    // Step 3: trim trailing SAM-D ops — these were WFA2-I trailing reference clips
    while( !cigar.empty() && (cigar.back() & 0xF) == CigarOp::D ) {
        ref_span -= static_cast<int>( cigar.back() >> 4 );
        cigar.pop_back();
    }

    int const ref_end = ref_begin + ref_span;

    // Step 4: sanity checks — catch any CIGAR/coordinate inconsistency
    if( query_span != qlen ) {
        throw std::logic_error(
            "wfa2: CIGAR query span " + std::to_string( query_span ) +
            " != query length "       + std::to_string( qlen )
        );
    }
    if( ref_end != cig->end_h ) {
        throw std::logic_error(
            "wfa2: ref_end " + std::to_string( ref_end ) +
            " != end_h "     + std::to_string( cig->end_h )
        );
    }

    return Wfa2Result{
        .ref_begin = ref_begin,
        .ref_end   = ref_end,
        .score     = cig->score,
        .status    = Wfa2Result::Status::Ok,
        .cigar     = std::move( cigar ),
    };
}

} // namespace spear::align
