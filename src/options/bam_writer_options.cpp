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

#include "options/bam_writer_options.hpp"

#include "tools/misc.hpp"

#include <stdexcept>

// =================================================================================================
//     Setup
// =================================================================================================

void BamWriterOptions::add_bam_opts_to_app(
    CLI::App* sub, std::string const& group
) {
    internal_check(
        !format.option && !compression.option && !multi_hit_mode.option && !max_hits.option &&
        !score_window.option && !keep_unmapped_opt.option && !read_group_id.option &&
        !read_group_sample.option && !read_group_library.option && !read_group_platform.option &&
        !read_group_raw.option,
        "Cannot use the same BamWriterOptions object multiple times."
    );

    // Output format.
    format = sub->add_option(
        "--bam-format",
        format.value,
        "Output file format: 'sam' (uncompressed text) or 'bam' (BGZF-compressed binary)."
    );
    format.option->transform( CLI::IsMember({ "sam", "bam" }, CLI::ignore_case ));
    format.option->group( group );

    // BGZF compression level.
    compression = sub->add_option(
        "--bam-compression",
        compression.value,
        "BGZF compression level, 0 (none) to 9 (max); ignored for --bam-format sam."
    );
    compression.option->check( CLI::Range( 0, 9 ));
    compression.option->group( group );

    // Multi-hit mode.
    multi_hit_mode = sub->add_option(
        "--bam-multi-hit-mode",
        multi_hit_mode.value,
        "How to write multiple alignment hits for a read: 'separate' (one BAM record per hit, "
        "non-primary hits flagged as secondary) or 'xa' (primary record only, with the remaining "
        "hits packed into an XA:Z: aux tag)."
    );
    multi_hit_mode.option->transform( CLI::IsMember({ "separate", "xa" }, CLI::ignore_case ));
    multi_hit_mode.option->group( group );

    // Max hits per read.
    max_hits = sub->add_option(
        "--bam-max-hits",
        max_hits.value,
        "Keep at most this many alignment hits per read (the strongest, by alignment score), "
        "discarding any beyond that count before writing. The default value of 0 means no limit."
    );
    max_hits.option->check( CLI::NonNegativeNumber );
    max_hits.option->group( group );

    // Score window.
    score_window = sub->add_option(
        "--bam-score-window",
        score_window.value,
        "Discard alignment hits whose score is more than this many units below the best hit for "
        "the same read. A negative value (the default) means no limit."
    );
    score_window.option->group( group );

    // Keep unmapped reads.
    keep_unmapped_opt = sub->add_flag(
        "--bam-keep-unmapped",
        keep_unmapped_opt.value,
        "Write an explicit unmapped record for reads with no surviving alignment, instead of "
        "omitting them from the output entirely. Off by default, since for the environmental/"
        "metagenomic samples this tool targets, the unmapped fraction against any single "
        "reference is often the large majority, not a small residual worth keeping by default."
    );
    keep_unmapped_opt.option->group( group );

    // Read group: individual fields.
    read_group_id = sub->add_option(
        "--bam-read-group-id",
        read_group_id.value,
        "Fixed read group id (RG ID) applied to every read regardless of which input file it "
        "came from. If not given, the id is instead auto-derived per input file from that "
        "file's base name, so that reads from different files end up in different read groups. "
        "Only meaningful if --bam-read-group-sample is also given."
    );
    read_group_id.option->group( group );

    read_group_sample = sub->add_option(
        "--bam-read-group-sample",
        read_group_sample.value,
        "Read group sample name (RG SM). Read groups are opt-in (matching bwa/bowtie2/minimap2): "
        "if this is not given, no read group is written at all, regardless of the other "
        "--bam-read-group-* options."
    );
    read_group_sample.option->group( group );

    read_group_library = sub->add_option(
        "--bam-read-group-library",
        read_group_library.value,
        "Read group library (RG LB). Optional; only meaningful if --bam-read-group-sample is "
        "also given."
    );
    read_group_library.option->group( group );

    read_group_platform = sub->add_option(
        "--bam-read-group-platform",
        read_group_platform.value,
        "Read group platform (RG PL). Only meaningful if --bam-read-group-sample is also given."
    );
    read_group_platform.option->group( group );

    // Read group: raw alternative, mutually exclusive with the individual fields above.
    read_group_raw = sub->add_option(
        "--bam-read-group-raw",
        read_group_raw.value,
        "Complete, pre-formatted read group header line, e.g. "
        "'@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA', for callers that already have a full read-group "
        "specification in hand (mirroring bwa's/minimap2's -R) rather than wanting it built "
        "field-by-field. Mutually exclusive with the other --bam-read-group-* options."
    );
    read_group_raw.option->group( group );
    read_group_raw.option->excludes( read_group_id.option );
    read_group_raw.option->excludes( read_group_sample.option );
    read_group_raw.option->excludes( read_group_library.option );
    read_group_raw.option->excludes( read_group_platform.option );
}

// =================================================================================================
//     Settings factory
// =================================================================================================

spear::align::BamWriterSettings BamWriterOptions::make_bam_writer_settings() const
{
    spear::align::BamWriterSettings settings;
    settings.format = ( format.value == "sam" )
        ? spear::align::BamWriterSettings::Format::Sam
        : spear::align::BamWriterSettings::Format::Bam;
    settings.bam_compression = compression.value;
    settings.multi_hit_mode  = ( multi_hit_mode.value == "xa" )
        ? spear::align::BamWriterSettings::MultiHitMode::XaTag
        : spear::align::BamWriterSettings::MultiHitMode::SeparateRecords;
    settings.hit_selection.max_hits     = max_hits.value;
    settings.hit_selection.score_window = score_window.value;
    return settings;
}

bool BamWriterOptions::keep_unmapped() const
{
    return keep_unmapped_opt.value;
}

// =================================================================================================
//     Read groups
// =================================================================================================

namespace {

// Replace literal two-character "\t" escapes with a real tab byte. SAM header lines are
// tab-delimited, and shells (esp. with single quotes) do not expand \t themselves, so without
// this, the common bwa/minimap2-style idiom `--bam-read-group-raw '@RG\tID:...\tSM:...'` would
// silently fail to parse as a valid header line. This is the only escape sequence unescaped
// here (not \n, \\, etc.): SAM header syntax has no use for any other escape.
std::string unescape_tabs( std::string const& s )
{
    std::string result;
    result.reserve( s.size() );
    for( std::size_t i = 0; i < s.size(); ++i ) {
        if( s[i] == '\\' && i + 1 < s.size() && s[i + 1] == 't' ) {
            result += '\t';
            ++i;
        } else {
            result += s[i];
        }
    }
    return result;
}

} // anonymous namespace

void BamWriterOptions::populate_read_groups(
    spear::align::BamHeader& header,
    std::vector<std::string> const& file_bases
) const {
    if( read_group_raw.was_set() ) {
        header.raw_read_groups.push_back( unescape_tabs( read_group_raw.value ));
        return;
    }
    if( !read_group_sample.was_set() ) {
        // Read groups are opt-in; nothing to do.
        return;
    }

    header.read_groups.reserve( read_group_id.was_set() ? 1 : file_bases.size() );
    if( read_group_id.was_set() ) {
        header.read_groups.push_back( spear::align::BamHeader::ReadGroup{
            read_group_id.value, read_group_sample.value,
            read_group_library.value, read_group_platform.value
        });
    } else {
        for( auto const& base : file_bases ) {
            header.read_groups.push_back( spear::align::BamHeader::ReadGroup{
                base, read_group_sample.value,
                read_group_library.value, read_group_platform.value
            });
        }
    }
}

std::optional<std::string_view> BamWriterOptions::read_group_id_for_file(
    std::size_t file_index,
    std::vector<std::string> const& file_bases
) const {
    if( read_group_raw.was_set() || !read_group_sample.was_set() ) {
        return std::nullopt;
    }
    if( read_group_id.was_set() ) {
        return std::string_view( read_group_id.value );
    }
    return std::string_view( file_bases.at( file_index ) );
}
