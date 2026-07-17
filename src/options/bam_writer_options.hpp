#ifndef SPEAR_OPTIONS_BAM_WRITER_H_
#define SPEAR_OPTIONS_BAM_WRITER_H_

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

#include "CLI/CLI.hpp"

#include "spear/align/bam_writer.hpp"

#include "tools/cli_option.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

// =================================================================================================
//     BamWriterOptions
// =================================================================================================

/**
 * @brief CLI binding for BamWriter configuration.
 *
 * Pure CLI wrapper: registers the BAM/SAM output options with a CLI11 subcommand and converts
 * them to a BamWriterSettings, plus a separate helper to populate a BamHeader's read group(s),
 * since that needs per-input-file information (base names) that this class does not itself hold.
 *
 * Embed one instance in a command's options struct and call add_bam_opts_to_app() once.
 */
class BamWriterOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    BamWriterOptions()  = default;
    ~BamWriterOptions() = default;

    BamWriterOptions( BamWriterOptions const& ) = default;
    BamWriterOptions( BamWriterOptions&& )      = default;

    BamWriterOptions& operator=( BamWriterOptions const& ) = default;
    BamWriterOptions& operator=( BamWriterOptions&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup
    // -------------------------------------------------------------------------

    /**
     * @brief Register all BAM/SAM output options with @p sub.
     *
     * Must not be called more than once on the same instance.
     */
    void add_bam_opts_to_app( CLI::App* sub, std::string const& group = "BAM Output" );

    // -------------------------------------------------------------------------
    //     Settings factory
    // -------------------------------------------------------------------------

    /**
     * @brief Build a BamWriterSettings (format, compression, multi-hit mode, hit selection)
     * from the current option values. Does not touch read groups; see populate_read_groups().
     */
    spear::align::BamWriterSettings make_bam_writer_settings() const;

    /**
     * @brief Whether unmapped reads should get an explicit BAM_FUNMAP record written, rather
     * than being silently absent from the output.
     */
    bool keep_unmapped() const;

    /**
     * @brief Populate @p header's read_groups or raw_read_groups according to the current
     * option values, using @p file_bases (one base name per input file, same order/indices as
     * the caller's own input file list) to derive per-file read group ids where applicable.
     *
     * Behavior:
     * - If --bam-read-group-raw was given: header.raw_read_groups = { that string }.
     * - Else if --bam-read-group-sample was not given: no read group is added at all (header is
     *   left untouched), matching the convention of bwa/bowtie2/minimap2: read groups are opt-in.
     * - Else: one ReadGroup per entry in @p file_bases, each with id = that file's base name,
     *   unless --bam-read-group-id was given, in which case a single ReadGroup with that fixed
     *   id is used for every file instead.
     *
     * @p header is expected to have empty read_groups/raw_read_groups on entry.
     */
    void populate_read_groups(
        spear::align::BamHeader& header,
        std::vector<std::string> const& file_bases
    ) const;

    /**
     * @brief The effective read group id that reads from input file @p file_index should be
     * tagged with (ReadRecord::tag_rg), or std::nullopt if no per-read RG tagging should happen:
     * either because read groups are not active (no --bam-read-group-sample), or because a raw
     * read-group line is being used instead (its own id is free-form text we do not parse back
     * out, so per-read RG tagging is not attempted in that case).
     *
     * The returned view points into either @p file_bases or this object's own read_group_id
     * option value; both must outlive the returned view, consistent with populate_read_groups().
     */
    std::optional<std::string_view> read_group_id_for_file(
        std::size_t file_index,
        std::vector<std::string> const& file_bases
    ) const;

    // -------------------------------------------------------------------------
    //     Option members
    // -------------------------------------------------------------------------

    // Output format: "sam" (uncompressed text) or "bam" (BGZF-compressed binary).
    CliOption<std::string> format = "bam";

    // BGZF compression level 0 (none) to 9 (max); ignored for SAM.
    CliOption<int> compression = 6;

    // How to write multiple hits for a read: "separate" (one record per hit, secondaries
    // flagged) or "xa" (primary record + secondaries packed into an XA:Z: aux tag).
    CliOption<std::string> multi_hit_mode = "separate";

    // Keep at most this many hits per read (after other filters); 0 = no limit.
    CliOption<int> max_hits = 0;

    // Discard hits whose AS is more than this many units below the best AS; negative = disabled.
    CliOption<int> score_window = -1;

    // If set, write an explicit BAM_FUNMAP record for reads with no surviving alignment, rather
    // than omitting them from the output entirely.
    CliOption<bool> keep_unmapped_opt = false;

    // Fixed read group id applied to every read regardless of source file; if not set (default),
    // the id is instead auto-derived per input file from that file's base name. Mutually
    // exclusive with --bam-read-group-raw.
    CliOption<std::string> read_group_id;

    // Read group sample name (SM). Presence of this option is what triggers read group output
    // at all (matching bwa/bowtie2/minimap2: read groups are opt-in, not attempted by default).
    // Mutually exclusive with --bam-read-group-raw.
    CliOption<std::string> read_group_sample;

    // Read group library (LB); optional. Mutually exclusive with --bam-read-group-raw.
    CliOption<std::string> read_group_library;

    // Read group platform (PL). Mutually exclusive with --bam-read-group-raw.
    CliOption<std::string> read_group_platform = "ILLUMINA";

    // Complete, pre-formatted "@RG\t..." line, for callers that already have a full read-group
    // specification in hand (mirroring bwa's/minimap2's -R) rather than wanting it built
    // field-by-field. Mutually exclusive with all of the read_group_* options above.
    CliOption<std::string> read_group_raw;

};

#endif // SPEAR_OPTIONS_BAM_WRITER_H_
