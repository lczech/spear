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

#include "commands/map/index.hpp"
#include "commands/map/format.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "spear/inverted_index/builder.hpp"
#include "spear/reference/collection.hpp"

#include "genesis/sequence/format/fasta_reader.hpp"
#include "genesis/sequence/format/fastx_input_stream.hpp"
#include "genesis/sequence/kmer/canonical_encoding.hpp"
#include "genesis/sequence/kmer/extractor.hpp"
#include "genesis/sequence/kmer/function.hpp"
#include "genesis/sequence/sequence.hpp"
#include "genesis/util/core/fs.hpp"
#include "genesis/util/core/info.hpp"
#include "genesis/util/core/logging.hpp"
#include "genesis/util/format/json/document.hpp"
#include "genesis/util/format/json/writer.hpp"
#include "genesis/util/io/input_source.hpp"
#include "genesis/util/text/string.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_map_index( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<MapIndexOptions>();
    auto sub = app.add_subcommand(
        "index",
        "Index one or more reference sequences for mapping short reads to them."
    );

    // -------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------

    // Input fasta files or directories.
    options->fasta.add_multi_file_input_opt_to_app(
        sub, "fasta", "fasta/fasta.gz",
        "(fasta|fas|fna|fa)(\\.gz)?", "{fasta,fas,fna,fa}[.gz]"
    );

    // -------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------

    // K-mer size.
    options->k = sub->add_option(
        "--k",
        options->k.value,
        "Size of the k-mers."
    );
    options->k.option->check( CLI::Range( size_t{1}, size_t{31} ));
    options->k.option->group( "Settings" );

    // Canonical k-mers.
    options->canonical = sub->add_flag(
        "--canonical",
        options->canonical.value,
        "Index canonical k-mers (the k-mer and its reverse complement), halving the number of "
        "distinct terms. If not set, k-mers are indexed as-is, and reverse-complement matches "
        "need to be found by separately querying the reverse complement k-mers of a read against "
        "the index, which can reduce alignment time."
    );
    options->canonical.option->group( "Settings" );

    // Genome bin width.
    options->genome_bin_width = sub->add_option(
        "--genome-bin-width",
        options->genome_bin_width.value,
        "Width of the bins (in bases) that genome positions are grouped into for the occurrences."
    );
    options->genome_bin_width.option->check( CLI::PositiveNumber );
    options->genome_bin_width.option->group( "Settings" );

    // -------------------------------------------------------------------
    //     Builder options
    // -------------------------------------------------------------------

    options->builder.add_build_opts_to_app( sub );

    // -------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------

    // Output options.
    options->output.set_group( "Output" );
    options->output.add_default_output_opts_to_app( sub );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( spear_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_map_index( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

namespace {

// -------------------------------------------------------------------------------------------------
//      run_map_index_impl
// -------------------------------------------------------------------------------------------------

/**
 * @brief Implementation of run_map_index(), templated on the position type used by the
 * inverted index, so that we can offer both 32 and 64 bit positions.
 */
template<typename PositionT>
void run_map_index_impl( MapIndexOptions const& options )
{
    using namespace genesis::sequence;
    using namespace genesis::util;

    // Check that none of our output files already exist, unless overwriting is allowed,
    // before running any expensive computations.
    options.output.check_output_files_nonexistence(
        std::string( "map-index" ), std::vector<std::string>{ "sidx", "json" }
    );

    // Prepare param shorthands.
    auto const k = static_cast<uint8_t>( options.k.value );
    auto const w = options.genome_bin_width.value;
    auto const canonical = options.canonical.value;

    // Determine the number of distinct terms (k-mers) that we index.
    std::size_t const num_term_indices = canonical
        ? number_of_canonical_kmers( k )
        : ( std::size_t{ 1 } << ( 2 * k ))
    ;

    // Set up the inverted index builder.
    using Builder = spear::inverted_index::InvertedIndexBuilder<PositionT>;
    Builder builder( num_term_indices, options.builder.make_builder_config<Builder>() );

    // Set up the canonical encoder, if requested. We use an optional-like raw pointer here,
    // as MinimalCanonicalEncoding has no default constructor.
    std::unique_ptr<MinimalCanonicalEncoding> encoder;
    if( canonical ) {
        encoder = std::make_unique<MinimalCanonicalEncoding>( k );
    }

    // Collect reference sequence metadata as we process each sequence.
    spear::reference::ReferenceCollection reference_col;

    // Get the thread pool to dispatch per-sequence work to.
    auto thread_pool = global_options.thread_pool();

    // -------------------------------------------------------------------
    //     Main loop
    // -------------------------------------------------------------------

    LOG_MSG2
        << "Memory usage before processing sequences: "
        << text::to_string_byte_format( core::info_process_current_memory_usage() )
    ;

    // Sequentially iterate over all input fasta files and their sequences, dispatching the k-mer
    // extraction and index building for each sequence to the thread pool. We keep track of the
    // total length of all sequences seen so far, which gives us the global offset of each new
    // sequence in the concatenated coordinate space that the index positions refer to.
    // Duplicate sequence names across files are allowed; they are disambiguated in the reference
    // catalog by their file index.
    std::size_t total_length = 0;
    std::size_t sequence_count = 0;
    auto const& fasta_paths = options.fasta.file_paths();
    for( auto const& fasta_path : fasta_paths ) {
        LOG_MSG1 << "Processing file " << fasta_path;

        // Absolute, symlink-resolved path recorded in the reference collection (and hence the
        // manifest), so that `spear map align` can later reliably reload this same file
        // regardless of its own working directory -- unlike fasta_path itself, which may be
        // relative to whatever directory this command happened to be run from.
        auto const canonical_fasta_path = core::real_path( fasta_path );

        auto fasta_input = FastaInputStream( io::from_file( fasta_path ));
        fasta_input.reader().validate_labels( false );
        for( auto& sequence : fasta_input ) {
            LOG_MSG2 << "Processing sequence " << sequence.label();
            ++sequence_count;
            auto const offset = total_length;
            auto const seq_length = sequence.sites().size();
            total_length += seq_length;

            // Check that the new total length still fits into the position type of the genome bins.
            auto constexpr max_pos = static_cast<std::size_t>( std::numeric_limits<PositionT>::max() );
            if( total_length > 0 && ( total_length - 1 ) / w > max_pos ) {
                throw std::runtime_error(
                    "Reference genome is too large for the configured genome bin width and position "
                    "type: the highest genome bin index (" + std::to_string(( total_length - 1 ) / w ) +
                    ") does not fit into " + std::to_string( sizeof(PositionT) * 8 ) + " bits. "
                    "Increase --genome-bin-width, or use a larger --position-bits value."
                );
            }

            // Record the sequence in the reference collection (length-only; sites not needed here).
            reference_col.add( canonical_fasta_path, sequence.label(), seq_length );

            // Process the k-mers of this sequence in a separate thread.
            auto seq_str = std::make_shared<std::string>( std::move( sequence.sites() ));
            thread_pool->enqueue_detached( [&builder, &encoder, k, w, offset, seq_str]()
            {
                auto extractor = KmerExtractor( k, std::move( *seq_str ));
                for( auto const& kmer : extractor ) {
                    auto const term_index = static_cast<typename Builder::term_index_type>(
                        encoder ? encoder->encode( kmer ) : kmer.value
                    );
                    auto const global_pos = offset + kmer.location;
                    auto const bin_index = static_cast<PositionT>( global_pos / w );
                    builder.add( term_index, bin_index );
                }
            });
        }
    }
    thread_pool->wait_for_all_pending_tasks();

    LOG_MSG2
        << "Memory usage after processing sequences: "
        << text::to_string_byte_format( core::info_process_current_memory_usage() )
    ;
    LOG_MSG
        << "Processed " << fasta_paths.size() << " file(s), "
        << sequence_count << " sequences, total length "
        << total_length << " bases."
    ;

    // -------------------------------------------------------------------
    //     Write output
    // -------------------------------------------------------------------

    // Make sure that the output directory exists before we write any files into it.
    genesis::util::core::dir_create( options.output.out_dir );

    // Write the inverted index itself. We do not need to call builder.finalize() here,
    // as write() already flushes all pending data as part of its single pass over the entries.
    auto const index_path = options.output.get_output_filename( "map-index", "sidx" );
    auto const stats = builder.write( index_path );
    LOG_MSG << "Wrote index to " << index_path;

    // Statistics output for the user.
    size_t posting_count_sum = 0;
    for( size_t i = 1; i + 1 < stats.posting_count_histogram.size(); ++i ) {
        posting_count_sum += stats.posting_count_histogram[i];
    }
    LOG_MSG1 << "Builder statistics:";
    LOG_MSG1 << "  Total blob size: " << text::to_string_byte_format( stats.total_blob_bytes );
    LOG_MSG1 << "  Empty k-mers:    " << stats.posting_count_histogram.front();
    LOG_MSG1 << "  Regular k-mers:  " << posting_count_sum;
    LOG_MSG1 << "  Capped k-mers:   " << stats.posting_count_histogram.back();
    LOG_MSG2 << "Detailed k-mer occurrence histogram (excluding empty and capped k-mers):";
    for( size_t i = 1; i + 1 < stats.posting_count_histogram.size(); ++i ) {
        if( stats.posting_count_histogram[i] > 0 ) {
            LOG_MSG2 << "  K-mers with " << i << " occurrences: " << stats.posting_count_histogram[i];
        }
    }

    // Write the JSON manifest with two sub-objects: "index" (build parameters) and
    // "reference" (sequence catalog for coordinate resolution at query time).
    MapIndexManifest manifest;
    manifest.metadata.sidx_path                = genesis::util::core::file_basename( index_path );
    manifest.metadata.k                        = k;
    manifest.metadata.canonical                = canonical;
    manifest.metadata.genome_bin_width         = w;
    manifest.metadata.max_occurrences_per_kmer = options.builder.max_postings_per_term.value;
    manifest.metadata.position_bits            = options.builder.position_bits.value;
    manifest.metadata.total_genome_length      = total_length;
    manifest.reference                         = std::move( reference_col );

    auto const json_target = options.output.get_output_target( "map-index", "json" );
    genesis::util::format::JsonWriter().write( manifest_to_json( manifest ), json_target );
}

} // anonymous namespace

void run_map_index( MapIndexOptions const& options )
{
    if( options.builder.position_bits.value == 64 ) {
        run_map_index_impl<std::uint64_t>( options );
    } else if( options.builder.position_bits.value == 32 ) {
        run_map_index_impl<std::uint32_t>( options );
    } else {
        throw std::invalid_argument(
            "Invalid --position-bits value " +
            std::to_string( options.builder.position_bits.value ) +
            "; must be 32 or 64."
        );
    }
}
