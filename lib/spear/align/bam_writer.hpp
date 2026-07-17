#ifndef SPEAR_LIB_ALIGN_BAM_WRITER_H_
#define SPEAR_LIB_ALIGN_BAM_WRITER_H_

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
#include <utility>
#include <vector>

#include <genesis/sequence/sequence_set.hpp>

#include "spear/align/read_record.hpp"

namespace spear::align {

// =================================================================================================
//     BamHeader
// =================================================================================================

/**
 * @brief SAM/BAM file header carrying reference sequences and optional metadata.
 *
 * Construct from a SequenceSet representing the full reference catalog. Names and lengths
 * are extracted and stored internally, so the SequenceSet need not remain alive after
 * construction. An optional read group can be added before passing to BamWriter.
 */
struct BamHeader
{
    /**
     * @brief SAM/BAM `@RG` (read group) header record.
     *
     * At minimum, @p id and @p sample must be set. @p library is optional.
     * Downstream tools such as mapDamage and GATK expect a read group to be present.
     */
    struct ReadGroup
    {
        std::string id;
        std::string sample;
        std::string library;
        std::string platform = "ILLUMINA";
    };

    /// Reference sequence names and lengths for @SQ header lines.
    std::vector<std::pair<std::string, int64_t>> references;

    /// Read groups written as @RG header lines, one per entry (e.g. one per input file, each
    /// with its own id). Mutually exclusive with @p raw_read_groups: BamWriter throws if both
    /// are non-empty.
    std::vector<ReadGroup> read_groups;

    /// Alternative to @p read_groups: complete, pre-formatted "@RG\t..." lines (one entry per
    /// line, as accepted by htslib's sam_hdr_add_lines()), for callers that already have a full
    /// read-group specification in hand -- e.g. mirroring bwa's/minimap2's `-R` -- rather than
    /// wanting it built field-by-field. Mutually exclusive with @p read_groups.
    std::vector<std::string> raw_read_groups;

    /**
     * @brief Default-construct an empty header; populate @p references manually
     * before passing to BamWriter.
     */
    BamHeader() = default;

    /**
     * @brief Construct from a SequenceSet; extracts label and length from each entry.
     */
    explicit BamHeader( genesis::sequence::SequenceSet const& refs );
};

// =================================================================================================
//     BamWriterSettings
// =================================================================================================

/**
 * @brief Configuration for a BamWriter instance.
 *
 * All settings are fixed at construction time. Defaults produce a standard compressed BAM
 * with one record per alignment hit and no hit filtering.
 */
struct BamWriterSettings
{
    /**
     * @brief Output file format.
     */
    enum class Format {
        Sam,   ///< Uncompressed SAM text.
        Bam,   ///< BGZF-compressed binary BAM.
        Cram,  ///< Reference-based CRAM (not yet implemented).
    };

    /**
     * @brief How secondary alignment hits in a ReadRecord are written.
     *
     * After HitSelectionPolicy filtering, hits[0] is always written as the primary record.
     * SeparateRecords writes each remaining hit as an additional BAM record with AlignmentFlags::Secondary.
     * XaTag packs all remaining hits into a single XA:Z: aux tag on the primary record.
     */
    enum class MultiHitMode {
        SeparateRecords, ///< One BAM record per hit; non-primary hits get AlignmentFlags::Secondary.
        XaTag,           ///< Primary record + secondaries packed into XA:Z: aux tag.
    };

    /// Output file format; default is BAM.
    Format       format          = Format::Bam;

    /// BGZF compression level 0 (none) to 9 (max); ignored for SAM and CRAM.
    int          bam_compression = 6;

    /// How to write multiple hits from a ReadRecord; default is one record per hit.
    MultiHitMode multi_hit_mode  = MultiHitMode::SeparateRecords;

    /// Controls which hits from a ReadRecord are written; see HitSelectionPolicy.
    HitSelectionPolicy hit_selection;
};

// =================================================================================================
//     BamWriter
// =================================================================================================

/**
 * @brief Writes ReadRecord objects to a SAM or BAM file.
 *
 * The output file is opened at construction and closed (and flushed) at destruction.
 * Not thread-safe. For ordered multi-threaded output, wrap write() as the output
 * function of a genesis::util::threading::SequentialOutputBuffer<ReadRecord>;
 * the caller assigns sequence IDs from a read counter when enqueuing to the thread pool.
 */
class BamWriter
{
public:

    /**
     * @brief Open @p path for writing with the given @p header and optional @p settings.
     *
     * @p header is taken by value and moved into the writer; pass `std::move(header)` to avoid
     * a copy when the header is no longer needed at the call site.
     * @throws std::runtime_error if the file cannot be opened or the header cannot be written.
     */
    BamWriter(
        std::string const& path,
        BamHeader header,
        BamWriterSettings settings = {}
    );
    ~BamWriter();

    BamWriter( BamWriter const& ) = delete;
    BamWriter& operator=( BamWriter const& ) = delete;
    BamWriter( BamWriter&& ) noexcept;
    BamWriter& operator=( BamWriter&& ) noexcept;

    /**
     * @brief Write all alignment hits from @p record according to the configured MultiHitMode.
     *
     * Hits are filtered by HitSelectionPolicy first. If no hits survive filtering, nothing
     * is written. hits[0] is treated as primary; the caller is responsible for ordering.
     * The record is consumed (moved from).
     */
    void write( ReadRecord&& record );

private:

    // Private implementation to keep htslib headers out of the public interface.
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace spear::align

#endif // SPEAR_LIB_ALIGN_BAM_WRITER_H_
