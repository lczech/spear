#ifndef SPEAR_LIB_ALIGN_WFA2_CIGAR_H_
#define SPEAR_LIB_ALIGN_WFA2_CIGAR_H_

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
#include <string_view>
#include <vector>

namespace spear::align {

// =================================================================================================
//     Wfa2CigarResult
// =================================================================================================

/**
 * @brief Result of parsing a WFA2 raw backtrace into a SAM-compliant CIGAR.
 */
struct Wfa2CigarResult
{
    /// 0-based inclusive alignment start in the reference window passed to WFA2.
    int ref_begin = 0;

    /// 0-based exclusive alignment end in the reference window passed to WFA2.
    int ref_end = 0;

    /// Edit distance (NM): mismatches + inserted/deleted bases.
    int edit_distance = 0;

    /// CIGAR in BAM uint32_t RLE format, as `(length << 4) | CigarOp::Op`.
    std::vector<uint32_t> cigar;
};

// =================================================================================================
//     build_sam_cigar_from_wfa2_ops
// =================================================================================================

/**
 * @brief Build a SAM-compliant CIGAR, ref_begin/ref_end, and edit distance directly from
 * WFA2's raw per-base backtrace operations.
 *
 * @p ops is WFA2's raw backtrace buffer (`cigar_t::operations`, sliced to
 * `[begin_offset, end_offset)`), using WFA2's pattern/text convention: 'M'/'X' are
 * diagonal (query- and reference-consuming) steps, 'I' is reference-only (SAM D), 'D' is
 * query-only (SAM I). WFA2 names operations from the pattern-matching viewpoint, which is
 * inverted relative to SAM's reference-to-query viewpoint. This convention follows
 * directly from the wavefront recursion itself (see WFA2's wavefront_backtrace.c), not
 * from a formatting choice in WFA2's SAM-CIGAR printing code, so it is expected to remain
 * stable even if that printing code changes upstream (https://github.com/smarco/WFA2-lib/pull/123).
 *
 * The aligner runs with free reference ends (fitting/semi-global alignment), so any
 * leading or trailing run of 'I' is reference padding outside the real alignment, not an
 * actual edit event; those runs are trimmed and folded into ref_begin/ref_end rather than
 * emitted as CIGAR operations.
 *
 * WFA2's 'M' only means "the match function accepted this step" -- under aDNA damage
 * tolerance, that includes tolerated C->T/G->A substitutions, not just literal sequence
 * equality. This function resolves that ambiguity itself by comparing @p query and
 * @p target directly, so the returned CIGAR and edit_distance are always correct
 * regardless of why WFA2 treated a given position as a match.
 *
 * @p target is the full reference window that was passed to the WFA2 aligner; positions
 * in @p ops index into it starting at h=0 (the very start of the window, including
 * through any leading padding).
 *
 * @p use_extended_cigar selects `=`/`X` (true) vs. plain `M` (false) for match/mismatch
 * operations; edit_distance is always computed regardless of this flag.
 */
Wfa2CigarResult build_sam_cigar_from_wfa2_ops(
    std::string_view ops,
    std::string_view query,
    std::string_view target,
    bool use_extended_cigar
);

} // namespace spear::align

#endif // SPEAR_LIB_ALIGN_WFA2_CIGAR_H_
