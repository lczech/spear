#ifndef SPEAR_LIB_ALIGN_CIGAR_H_
#define SPEAR_LIB_ALIGN_CIGAR_H_

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

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace spear::align {

// =================================================================================================
//     CIGAR utilities
// =================================================================================================

/**
 * @brief Convert a CIGAR vector to its SAM text form (e.g. "36M2I5D").
 *
 * Each uint32_t element encodes one run as `(length << 4) | op_code`, where op codes are:
 *
 * ```
 * 0=M  1=I  2=D  3=N  4=S  5=H  6=P  7==  8=X
 * ```
 *
 * This is the binary BAM CIGAR format, also used by WFA2's cigar_t::cigar_buffer.
 */
inline std::string cigar_to_string( std::vector<uint32_t> const& cigar )
{
    static constexpr char ops[] = "MIDNSHP=X";
    std::string result;
    result.reserve( cigar.size() * 4 );
    for( uint32_t const c : cigar ) {
        result += std::to_string( c >> 4 );
        result += ops[ c & 0xF ];
    }
    return result;
}

/**
 * @brief Compute the edit distance (NM tag) from an extended (X/=) CIGAR vector.
 *
 * Edit distance = sum of X (mismatch), I (insertion), and D (deletion) op lengths.
 * Only valid for extended CIGAR where `=` encodes a sequence match and `X` a mismatch.
 *
 * @throws std::invalid_argument if any `M`(0) op is encountered; use the two-sequence
 *         overload for M-CIGAR, which requires the read and reference sequences.
 */
inline int edit_distance_from_cigar( std::vector<uint32_t> const& cigar )
{
    // BAM op codes: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
    int dist = 0;
    for( uint32_t const c : cigar ) {
        uint32_t const op  = c & 0xF;
        uint32_t const len = c >> 4;
        if( op == 7 ) { // =: sequence match, no contribution
            (void) len;
        } else if( op == 8 || op == 1 || op == 2 ) { // X, I, D: contribute to edit distance
            dist += static_cast<int>( len );
        } else if( op == 0 ) { // M: ambiguous, cannot compute edit distance without sequences
            throw std::invalid_argument(
                "edit_distance_from_cigar: M op requires sequence context; "
                "use the two-sequence overload or set use_extended_cigar = true"
            );
        } else { // N(3), S(4), H(5), P(6): unexpected in gap-affine alignment output
            throw std::invalid_argument(
                "edit_distance_from_cigar: unexpected CIGAR op code " + std::to_string( op )
            );
        }
    }
    return dist;
}

/**
 * @brief Compute the edit distance (NM tag) from a CIGAR vector, handling both M and X/= formats.
 *
 * Processes each op as encountered:
 * - `=`(7): sequence match — no contribution.
 * - `X`(8): sequence mismatch — contributes length to edit distance.
 * - `I`(1): insertion in query — contributes length.
 * - `D`(2): deletion from query — contributes length.
 * - `M`(0): ambiguous — compares @p query and @p ref_window character by character
 *            to count mismatches; then adds insertion/deletion lengths from I/D ops.
 *
 * @p ref_window is the reference sequence that was passed to the aligner (the full window).
 * @p ref_begin is the 0-based start of the alignment within that window (from the alignment result).
 */
inline int edit_distance_from_cigar(
    std::vector<uint32_t> const& cigar,
    std::string_view              query,
    std::string_view              ref_window,
    int                           ref_begin
) {
    // BAM op codes: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
    int dist  = 0;
    int q_pos = 0;
    int r_pos = ref_begin;
    for( uint32_t const c : cigar ) {
        uint32_t const op  = c & 0xF;
        int      const len = static_cast<int>( c >> 4 );
        if( op == 7 ) { // =: sequence match
            q_pos += len;
            r_pos += len;
        } else if( op == 8 ) { // X: sequence mismatch
            dist  += len;
            q_pos += len;
            r_pos += len;
        } else if( op == 1 ) { // I: insertion in query
            dist  += len;
            q_pos += len;
        } else if( op == 2 ) { // D: deletion from query
            dist  += len;
            r_pos += len;
        } else if( op == 0 ) { // M: compare sequences to resolve matches vs mismatches
            for( int i = 0; i < len; ++i ) {
                if( query[q_pos] != ref_window[r_pos] ) {
                    ++dist;
                }
                ++q_pos;
                ++r_pos;
            }
        } else { // N(3), S(4), H(5), P(6): unexpected in gap-affine alignment output
            throw std::invalid_argument(
                "edit_distance_from_cigar: unexpected CIGAR op code " + std::to_string( op )
            );
        }
    }
    return dist;
}

} // namespace spear::align

#endif // SPEAR_LIB_ALIGN_CIGAR_H_
