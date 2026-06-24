#ifndef SPEAR_LIB_ALIGN_READ_RECORD_H_
#define SPEAR_LIB_ALIGN_READ_RECORD_H_

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
#include <optional>
#include <string>
#include <vector>

#include <genesis/sequence/sequence.hpp>

namespace spear::align {

// =================================================================================================
//     AlignmentFlags
// =================================================================================================

/**
 * @brief Bitmask constants for AlignmentHit::flags.
 *
 * Values match the SAM/BAM FLAG field specification; the underlying concepts (strand, pairing,
 * secondary) are universal to most alignment formats.
 */
namespace AlignmentFlags {
    // Static assertions in bam_writer.cpp verify at compile time
    // that these values agree with the htslib BAM_F* constants.
    enum AlignmentFlag : uint16_t {
        Paired        = 0x001, ///< Read is paired in sequencing.
        ProperPair    = 0x002, ///< Read and mate align in expected orientation and distance.
        Unmapped      = 0x004, ///< Read is unmapped.
        MateUnmapped  = 0x008, ///< Mate is unmapped.
        Reverse       = 0x010, ///< Read aligns to the reverse strand.
        MateReverse   = 0x020, ///< Mate aligns to the reverse strand.
        Read1         = 0x040, ///< First read in a pair.
        Read2         = 0x080, ///< Second read in a pair.
        Secondary     = 0x100, ///< Secondary alignment.
        QcFail        = 0x200, ///< Read failed platform/vendor quality checks.
        Duplicate     = 0x400, ///< PCR or optical duplicate.
        Supplementary = 0x800, ///< Supplementary (chimeric) alignment.
    };
} // namespace AlignmentFlags

// =================================================================================================
//     HitSelectionPolicy
// =================================================================================================

/**
 * @brief Policy controlling which alignment hits from a ReadRecord are written to output.
 *
 * All three filters are independent and applied in combination: a hit must pass every active
 * filter to be written. A zero or negative value disables the respective filter.
 *
 * The score_window filter operates on AlignmentHit::tag_as (alignment score): a hit is retained
 * if its score is within @p score_window penalty units of the best score among all hits
 * in the record. For example, a score_window of 5 retains hits with AS ≥ (best_AS - 5).
 */
struct HitSelectionPolicy
{
    /// Keep at most this many hits (after other filters); 0 = no limit.
    int max_hits     = 0;

    /// Discard hits with MAPQ below this threshold; 0 = no filter.
    int min_mapq     = 0;

    /// Discard hits whose AS is more than this many units below the best AS; negative = no filter.
    int score_window = -1;
};

// =================================================================================================
//     AlignmentHit
// =================================================================================================

/**
 * @brief One alignment of a read to a specific reference location.
 *
 * Owns all data so it can be moved safely across threads.
 * For single-end reads leave the mate_* fields at their defaults and do not set
 * AlignmentFlags::Paired in @p flags. For paired-end reads set AlignmentFlags::Paired,
 * AlignmentFlags::Read1 or AlignmentFlags::Read2, and fill the mate fields accordingly.
 */
struct AlignmentHit
{
    // Fields are ordered for minimal struct padding (120 bytes on 64-bit): core alignment fields
    // first, then paired-end mate fields, then quality annotation tags.

    // -------------------------------------------------------------------------
    //     Core: where and how the read maps
    // -------------------------------------------------------------------------

    /// 0-based leftmost mapped position on the reference.
    int64_t               ref_pos      = 0;

    /// 0-based index into BamHeader::references; -1 = unmapped.
    int32_t               ref_id       = -1;

    /// Alignment flags; use AlignmentFlags::AlignmentFlag values.
    uint16_t              flags        = 0;

    /// Mapping quality (255 = not available).
    uint8_t               mapq         = 255;

    /// CIGAR in BAM uint32_t RLE format (`(length << 4) | op_code`); empty for unmapped reads.
    std::vector<uint32_t> cigar;

    // -------------------------------------------------------------------------
    //     Mate: paired-end fields
    // -------------------------------------------------------------------------

    /// Mate 0-based leftmost position; 0 if no mate.
    int64_t               mate_ref_pos = 0;

    /// Observed template length (insert size); 0 if not applicable.
    int64_t               template_len = 0;

    /// Mate chromosome index; -1 if no mate or mate is unmapped.
    int32_t               mate_ref_id  = -1;

    // -------------------------------------------------------------------------
    //     Tags: alignment quality annotations
    // -------------------------------------------------------------------------

    /// AS aux tag: alignment score from the aligner (always written).
    int                   tag_as       = 0;

    /// NM aux tag: edit distance (mismatches + indel bases); omitted if empty.
    std::optional<int>    tag_nm;

    /// MD aux tag: mismatch/deletion reference string; omitted if empty.
    /// Requires the reference window at computation time; use samtools calmd if not set.
    std::optional<std::string> tag_md;
};

// =================================================================================================
//     ReadRecord
// =================================================================================================

/**
 * @brief A sequencing read together with all of its alignment hits.
 *
 * Owns all data so it can be moved safely across threads.
 *
 * @p hits[0] is treated as the primary alignment by output writers; the caller is responsible
 * for ordering hits with the preferred primary first. Output format and secondary-hit handling
 * (e.g., separate BAM records vs XA tag) are decided by the writer at serialisation time.
 */
struct ReadRecord
{
    /// Read name, nucleotide sequence, and optional Phred quality scores (empty = unavailable).
    genesis::sequence::Sequence read;

    /// RG aux tag: read group ID; should match BamHeader::read_group::id when set.
    std::optional<std::string> tag_rg;

    /// All alignment hits for this read; [0] is primary.
    std::vector<AlignmentHit> hits;
};

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

} // namespace spear::align

#endif // SPEAR_LIB_ALIGN_READ_RECORD_H_
