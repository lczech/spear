#ifndef SPEAR_LIB_ALIGN_WFA2_H_
#define SPEAR_LIB_ALIGN_WFA2_H_

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
#include <memory>
#include <string>
#include <vector>

namespace spear::align {

// =================================================================================================
//     Settings
// =================================================================================================

/**
 * @brief Configuration for a Wfa2Aligner instance.
 *
 * All settings are fixed at construction time. The aligner performs semi-global
 * (fitting) alignment: the query (read) is fully consumed, while both ends of the
 * target (reference window) are free. Gap-affine penalties follow the
 * Smith-Waterman-Gotoh model with a fixed match cost of 0.
 *
 * aDNA damage tolerance is enabled by setting non-zero zone sizes: bases within
 * @p damage_ct_end of the 5' end and @p damage_ga_start of the 3' end use a
 * custom match function that treats C→T and G→A substitutions as free matches,
 * respectively. When both zone sizes are zero the standard match function is used.
 */
struct Wfa2Settings
{
    /**
     * @brief Scope of alignment results to compute.
     *
     * Score computes the alignment score and ref_end only (no backtrace, faster).
     * Cigar additionally computes ref_begin and the full CIGAR string.
     */
    enum class Scope {
        Score, ///< `score` + `ref_end` only; faster, no backtrace.
        Cigar, ///< `score` + `ref_begin` + `ref_end` + CIGAR.
    };

    /// Scope of results to compute; fixed at construction.
    Scope scope = Scope::Cigar;

    /// Mismatch penalty (non-negative; match is always 0).
    int mismatch   = 4;
    /// Gap-open penalty (non-negative).
    int gap_open   = 4;
    /// Gap-extension penalty (non-negative).
    int gap_extend = 2;

    /// Number of bases from the 5' end where C→T is treated as a free match.
    /// Set to 0 to disable (default).
    int damage_ct_end   = 0;
    /// Number of bases from the 3' end where G→A is treated as a free match.
    /// Set to 0 to disable (default).
    int damage_ga_start = 0;

    /**
     * @brief Maximum alignment score (WFA steps) before aborting.
     *
     * Acts as a safety cap against pathological inputs. For reads ≤ 256 bp
     * the worst-case score is 256 × mismatch = 1024, so the default of 4096
     * guarantees completion at any realistic divergence.
     */
    int max_steps = 4096;

    /**
     * @brief Emit `X/=` CIGAR operators instead of the ambiguous `M` operator.
     *
     * When true (default), `=` marks a sequence match and `X` marks a mismatch.
     * When false, both are collapsed to `M` for compatibility with old tools.
     *
     * X/= CIGAR is preferred: it lets NM be computed from the CIGAR alone and
     * gives mapDamage and similar tools explicit mismatch positions. M-CIGAR
     * requires the reference window to derive NM (via samtools calmd).
     */
    bool use_extended_cigar = true;
};

// =================================================================================================
//     Result
// =================================================================================================

/**
 * @brief Result of a single WFA2 alignment call.
 */
struct Wfa2Result
{
    /**
     * @brief Alignment outcome.
     */
    enum class Status {
        Ok,               ///< Alignment completed successfully.
        ScoreCapExceeded, ///< Alignment score exceeded Wfa2Settings::max_steps.
        OutOfMemory,      ///< WFA2 hit its memory limit.
        Failed,           ///< Any other failure (unreachable end, etc.).
    };

    /// 0-based inclusive start in reference. Only valid for Scope::Cigar.
    int ref_begin = 0;

    /// 0-based exclusive end in reference.
    int ref_end = 0;

    /// WFA2 penalty score (≤0; 0 = perfect match, more negative = worse).
    int score = 0;

    /// Alignment success status.
    Status status = Status::Ok;

    /**
     * @brief CIGAR in BAM uint32_t RLE format, owned by this result.
     *
     * Only populated for Scope::Cigar and Status::Ok; empty otherwise.
     * Safe to move into AlignmentHit::cigar without any copy.
     */
    std::vector<uint32_t> cigar;
};

// =================================================================================================
//     Aligner
// =================================================================================================

/**
 * @brief WFA2 semi-global aligner (pattern fully consumed, text ends free).
 *
 * Not thread-safe. Create one instance per thread.
 * Scope (Score vs Cigar) and all settings are fixed at construction.
 */
class Wfa2Aligner
{
public:

    explicit Wfa2Aligner( Wfa2Settings settings );
    ~Wfa2Aligner();

    Wfa2Aligner( Wfa2Aligner const& ) = delete;
    Wfa2Aligner& operator=( Wfa2Aligner const& ) = delete;
    Wfa2Aligner( Wfa2Aligner&& ) noexcept;
    Wfa2Aligner& operator=( Wfa2Aligner&& ) noexcept;

    /// Run score-only alignment; returns score and ref_end, but not ref_begin or CIGAR.
    /// query is the sequencing read (pattern, fully consumed);
    /// target is the reference genome window (text, free ends).
    /// Throws if constructed with Scope::Cigar.
    Wfa2Result align_score( std::string const& query, std::string const& target );

    /// Run full alignment with backtrace; returns score, ref_begin, ref_end, and CIGAR.
    /// query is the sequencing read (pattern, fully consumed);
    /// target is the reference genome window (text, free ends).
    /// Throws if constructed with Scope::Score.
    Wfa2Result align_cigar( std::string const& query, std::string const& target );

private:

    // Private implementation idiom to hide WFA2's C API and internal state.
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace spear::align

#endif // SPEAR_LIB_ALIGN_WFA2_H_
