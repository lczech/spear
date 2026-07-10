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

#include "spear/align/wfa2_cigar.hpp"
#include "spear/align/cigar.hpp"

#include <cassert>
#include <cstddef>

namespace spear::align {

Wfa2CigarResult build_sam_cigar_from_wfa2_ops(
    std::string_view ops,
    std::string_view query,
    std::string_view target,
    bool use_extended_cigar
) {
    Wfa2CigarResult result;

    // Append one base of the given SAM op code, merging into the previous run if it matches.
    auto append_base_ = []( std::vector<uint32_t>& cigar, uint32_t op )
    {
        if( !cigar.empty() && (cigar.back() & 0xF) == op ) {
            cigar.back() += (uint32_t{ 1 } << 4);
        } else {
            cigar.push_back( (uint32_t{ 1 } << 4) | op );
        }
    };

    // Trim leading and trailing reference-only padding (free-end alignment clips),
    // symmetrically, directly on the raw WFA2 op chars, before any I/D swap.
    // Both are guaranteed to be a single contiguous run: pattern_begin_free=0 and
    // pattern_end_free=0 force full query consumption, so the optimal alignment never
    // pays for a real query-consuming op where free reference padding is available.
    std::size_t begin = 0;
    std::size_t end   = ops.size();
    while( begin < end && ops[begin] == 'I' ) {
        ++result.ref_begin;
        ++begin;
    }
    while( end > begin && ops[end - 1] == 'I' ) {
        --end;
    }

    // Single forward pass over the real alignment: swap WFA2's pattern/text convention
    // to SAM (I<->D), resolve 'M' into '='/X by comparing sequences, and RLE-encode.
    int q = 0;
    int r = result.ref_begin;
    for( std::size_t i = begin; i < end; ++i ) {
        char const raw = ops[i];
        if( raw == 'M' || raw == 'X' ) {
            assert( static_cast<std::size_t>( q ) < query.size() );
            assert( static_cast<std::size_t>( r ) < target.size() );
            bool const is_match = ( raw == 'M' ) && ( query[q] == target[r] );
            if( !is_match ) {
                ++result.edit_distance;
            }
            uint32_t const op = use_extended_cigar
                ? ( is_match ? CigarOp::E : CigarOp::X )
                : CigarOp::M;
            append_base_( result.cigar, op );
            ++q;
            ++r;
        } else if( raw == 'I' ) {
            // WFA2 raw I: reference-consuming only -> SAM D
            ++result.edit_distance;
            append_base_( result.cigar, CigarOp::D );
            ++r;
        } else if( raw == 'D' ) {
            // WFA2 raw D: query-consuming only -> SAM I
            ++result.edit_distance;
            append_base_( result.cigar, CigarOp::I );
            ++q;
        } else {
            // WFA2 gap-affine backtrace only ever emits M/X/I/D; anything else means our
            // understanding of WFA2's contract has changed and needs re-checking.
            assert( false && "build_sam_cigar_from_wfa2_ops: unexpected WFA2 op code" );
        }
    }

    result.ref_end = r;
    assert(
        q == static_cast<int>( query.size() ) &&
        "build_sam_cigar_from_wfa2_ops: CIGAR query span does not match query length"
    );

    return result;
}

} // namespace spear::align
