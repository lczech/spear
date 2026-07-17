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

#include "commands/map/align.hpp"
#include "commands/map/format.hpp"
#include "options/global.hpp"
#include "tools/cli_setup.hpp"

#include "spear/align/bam_writer.hpp"
#include "spear/align/cigar.hpp"
#include "spear/align/read_record.hpp"
#include "spear/align/wfa2.hpp"
#include "spear/inverted_index/index.hpp"
#include "spear/reference/collection.hpp"
#include "spear/reference/functions.hpp"
#include "spear/seeding/seed_filter.hpp"

#include "genesis/sequence/format/fastx_input_view_stream.hpp"
#include "genesis/sequence/function/code.hpp"
#include "genesis/sequence/kmer/canonical_encoding.hpp"
#include "genesis/sequence/kmer/extractor.hpp"
#include "genesis/util/core/logging.hpp"
#include "genesis/util/format/json/document.hpp"
#include "genesis/util/format/json/reader.hpp"
#include "genesis/util/io/input_source.hpp"
#include "genesis/util/threading/sequential_output_buffer.hpp"

#include <atomic>
#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// =================================================================================================
//      Setup
// =================================================================================================

void setup_map_align( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<MapAlignOptions>();
    auto sub = app.add_subcommand(
        "align",
        "Align (aka map) short reads to a reference genome."
    );

    // -------------------------------------------------------------------
    //     Input
    // -------------------------------------------------------------------

    // Index file and load mode.
    options->index_reader.add_reader_opts_to_app( sub, "Input" );

    // Input reads.
    options->reads.add_multi_file_input_opt_to_app(
        sub, "reads", "fasta/fastq",
        "(fa|fasta|fq|fastq)(\\.gz)?", "{fa,fasta,fq,fastq}[.gz]",
        true, "Input"
    );

    // -------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------

    options->kmer_seeding.add_kmer_seeding_opts_to_app( sub );
    options->wfa2.add_wfa2_opts_to_app( sub );
    options->bam_writer.add_bam_opts_to_app( sub );

    // -------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------

    // Output options.
    options->output.set_group( "Output" );
    options->output.add_default_output_opts_to_app( sub );

    options->emit_tsv = sub->add_flag(
        "--emit-tsv",
        options->emit_tsv.value,
        "Additionally write a diagnostic TSV alongside the BAM/SAM output, with seeding-stage "
        "information (e.g. number of candidate intervals found/kept, strand tie-break flags) "
        "that has no natural home in BAM. Off by default; BAM/SAM is the primary output."
    );
    options->emit_tsv.option->group( "Output" );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub->callback( spear_cli_callback(
        sub,
        {
            // Citation keys as needed
        },
        [ options ]() {
            run_map_align( *options );
        }
    ));
}

// =================================================================================================
//      Run
// =================================================================================================

namespace {


// -------------------------------------------------------------------------------------------------
//      run_map_align_impl
// -------------------------------------------------------------------------------------------------

/**
 * @brief Implementation of run_map_align(), templated on the position type used by the
 * inverted index, so that we can read both 32 and 64 bit indices.
 */
template<typename PositionT>
void run_map_align_impl(
    MapAlignOptions const& options,
    MapIndexMetadata const& meta,
    spear::reference::ReferenceCollection const& reference,
    std::string const& sidx_path
) {
    using namespace genesis::sequence;
    using namespace genesis::util;

    auto const k         = meta.k;
    auto const canonical = meta.canonical;
    auto const w         = meta.genome_bin_width;

    if( reference.empty() ) {
        throw std::runtime_error(
            "Index '" + options.index_reader.index_file.value + "' has no reference sequences."
        );
    }

    // Open the inverted index.
    using Index = spear::inverted_index::InvertedIndex<PositionT>;
    LOG_MSG << "Loading index " << sidx_path << " with genome bin width " << w;
    auto const mode = options.index_reader.is_pread()
        ? Index::OpenMode::kPread
        : Index::OpenMode::kLoadAll;
    Index index;
    index.open( sidx_path, mode );
    LOG_MSG << "Loaded index with " << index.num_terms() << " terms.";
    global_options.report_memory_usage();

    // Construct the seeding engine (validates bin_width, holds index reference).
    auto seeding = options.kmer_seeding.make_engine( index, w );

    // Build the seed filter config once; captured by reference into each per-read task below,
    // so that each of the many enqueued closures only stores a pointer, not a copy.
    auto const seed_filter_cfg = options.kmer_seeding.make_seed_filter_config();

    // Same for the WFA2 alignment settings; each thread builds its own Wfa2Aligner from this
    // (Wfa2Aligner is not thread-safe, so it cannot be shared across threads like seeding is).
    auto const wfa2_settings = options.wfa2.make_settings();

    // Set up the canonical encoder, if requested. We use a unique_ptr here,
    // as MinimalCanonicalEncoding has no default constructor.
    std::unique_ptr<MinimalCanonicalEncoding> encoder;
    if( canonical ) {
        encoder = std::make_unique<MinimalCanonicalEncoding>( k );
    }

    // -------------------------------------------------------------------
    //     Output setup
    // -------------------------------------------------------------------

    auto const bam_ext = ( options.bam_writer.format.value == "sam" )
        ? std::string( "sam" ) : std::string( "bam" );
    if( options.emit_tsv.value ) {
        options.output.check_output_files_nonexistence(
            std::vector<std::pair<std::string, std::string>>{
                { "map-align", bam_ext }, { "map-align", "tsv" }
            }
        );
    } else {
        options.output.check_output_files_nonexistence(
            std::string( "map-align" ), bam_ext
        );
    }
    genesis::util::core::dir_create( options.output.out_dir );

    // Per-input-file base names, computed once and kept stable for the whole run: ReadRecord::
    // tag_rg (below) is a non-owning string_view into these entries, so this vector must outlive
    // every buffered ReadRecord (it does -- function scope, past the thread_pool wait below), and
    // must not be resized/reallocated after this point (it never is).
    auto const file_bases = options.reads.base_file_names();

    // Build the BAM/SAM header: reference sequences, plus read group(s) if the user opted in
    // (populate_read_groups() is a no-op if --bam-read-group-sample was not given).
    spear::align::BamHeader bam_header;
    bam_header.references.reserve( reference.size() );
    for( std::size_t i = 0; i < reference.size(); ++i ) {
        bam_header.references.emplace_back(
            reference[i].name(), static_cast<int64_t>( reference[i].length )
        );
    }
    options.bam_writer.populate_read_groups( bam_header, file_bases );

    spear::align::BamWriter bam_writer(
        options.output.get_output_filename( "map-align", bam_ext ),
        std::move( bam_header ),
        options.bam_writer.make_bam_writer_settings()
    );
    bool const keep_unmapped = options.bam_writer.keep_unmapped();

    // Buffer that re-establishes the original (per-file, sequential) read order in the output,
    // even though reads are processed out of order across the thread pool.
    threading::SequentialOutputBuffer<spear::align::ReadRecord> bam_out_buff(
        [&]( spear::align::ReadRecord&& rr ) {
            bam_writer.write( std::move( rr ));
        }
    );

    // Optional diagnostic TSV alongside the BAM/SAM output; off by default (MapAlignOptions::
    // emit_tsv). Built via a unique_ptr so it is only constructed (and only ever written to)
    // when actually requested.
    std::shared_ptr<genesis::util::io::BaseOutputTarget> tsv_target;
    std::unique_ptr<threading::SequentialOutputBuffer<std::string>> tsv_out_buff;
    bool const emit_tsv = options.emit_tsv.value;
    if( emit_tsv ) {
        tsv_target = options.output.get_output_target( "map-align", "tsv" );
        tsv_target->ostream()
            << "read_file\tread_name\tread_length\tnum_kmers\tnum_seeds_found\t"
            << "num_seeds_kept\tcontig\tstart\tend\tpeak_hits\tstrand\talign_score\t"
            << "align_score_other_strand\tstrand_tie\tcigar\n"
        ;
        tsv_out_buff = std::make_unique<threading::SequentialOutputBuffer<std::string>>(
            [&]( std::string&& line ) {
                tsv_target->ostream() << line;
            }
        );
    }
    global_options.report_memory_usage();

    // -------------------------------------------------------------------
    //     Main loop
    // -------------------------------------------------------------------

    auto thread_pool = global_options.thread_pool();
    std::size_t total_reads = 0;
    auto const& read_files = options.reads.file_paths();

    LOG_MSG << "Processing " << read_files.size() << " read file(s)";
    for( std::size_t fi = 0; fi < read_files.size(); ++fi ) {
        auto const& path = read_files[fi];
        // Non-owning view into the stable file_bases vector built above; safe to copy by value
        // into per-read lambdas below, since the vector itself outlives them (see its own comment).
        std::string_view const file_base = file_bases[fi];
        auto const rg_id = options.bam_writer.read_group_id_for_file( fi, file_bases );

        for( auto const& rec : FastxInputViewStream( io::from_file( path ))) {
            // Extract the read name (first whitespace-delimited token of the label),
            // and copy the sequence, as the views are invalidated on the next increment.
            auto const& label = rec.label();
            auto const name_end = label.find_first_of( " \t" );
            std::string read_name(
                name_end == std::string_view::npos ? label : label.substr( 0, name_end )
            );
            auto seq_ptr = std::make_shared<std::string>( rec.sites() );

            // The index of this read in the original input order,
            // used to re-establish that order in the output.
            auto const read_idx = total_reads;
            ++total_reads;

            thread_pool->enqueue_detached(
                [&reference, &encoder, &bam_out_buff, &tsv_out_buff, &seeding, &seed_filter_cfg,
                 &wfa2_settings, canonical, keep_unmapped, emit_tsv, k, w, read_idx, file_base,
                 rg_id, read_name, seq_ptr]()
            {
                auto& seq = *seq_ptr;
                auto const read_length = seq.size();
                std::size_t const num_kmers = ( read_length >= k ) ? ( read_length - k + 1 ) : 0;

                std::ostringstream line;

                spear::align::ReadRecord rr;
                rr.read = genesis::sequence::Sequence( read_name, seq );
                if( rg_id ) {
                    rr.tag_rg = *rg_id;
                }

                using KmerSeedingT = spear::seeding::KmerSeeding<PositionT>;
                using SeedInterval = typename KmerSeedingT::SeedInterval;
                thread_local std::vector<SeedInterval> seed_intervals;
                std::size_t num_seeds_found = 0;
                std::size_t num_seeds_kept  = 0;

                if( canonical ) {
                    // Canonical index: a single seeding query already matches k-mers regardless
                    // of the read's orientation (MinimalCanonicalEncoding maps a k-mer and its
                    // reverse complement to the same term), so which strand actually matched is
                    // resolved later, at the alignment step below.
                    auto query = seeding.start_query();
                    KmerExtractor extractor( k, std::string_view( seq ));
                    for( auto const& kmer : extractor ) {
                        auto const term_index = static_cast<typename Index::term_index_type>(
                            encoder->encode( kmer )
                        );
                        if( !query.add_kmer( term_index ) ) {
                            break;  // read is truncated; remaining k-mers are not seeded
                        }
                    }
                    query.finish_query( read_length, seed_intervals );
                    num_seeds_found = seed_intervals.size();
                } else {
                    // Non-canonical index: terms are raw (non-canonical) k-mer values, so a
                    // reverse-strand read only matches reference terms via its reverse
                    // complement's k-mers, requiring a second seeding query. A single
                    // KmerExtractor pass over the forward read already yields both a k-mer's
                    // value and its rev_comp incrementally at no extra cost (see
                    // genesis::sequence::Kmer::rev_comp), so no second extraction pass over a
                    // materialized reverse-complement string is needed here -- only two
                    // downstream seeding queries, run sequentially (KmerSeeding allows only one
                    // active Query per thread at a time).
                    thread_local std::vector<typename Index::term_index_type> fwd_kmer_indices;
                    thread_local std::vector<typename Index::term_index_type> rev_kmer_indices;
                    fwd_kmer_indices.clear();
                    rev_kmer_indices.clear();

                    KmerExtractor extractor( k, std::string_view( seq ));
                    for( auto const& kmer : extractor ) {
                        if( fwd_kmer_indices.size() >= KmerSeedingT::kMaxKmersPerRead ) {
                            break;  // read is truncated; remaining k-mers are not seeded
                        }
                        fwd_kmer_indices.push_back(
                            static_cast<typename Index::term_index_type>( kmer.value )
                        );
                        rev_kmer_indices.push_back(
                            static_cast<typename Index::term_index_type>( kmer.rev_comp )
                        );
                    }

                    thread_local std::vector<SeedInterval> fwd_intervals;
                    thread_local std::vector<SeedInterval> rev_intervals;
                    seeding.run_query( fwd_kmer_indices, read_length, fwd_intervals );
                    seeding.run_query( rev_kmer_indices, read_length, rev_intervals );
                    for( auto& iv : fwd_intervals ) { iv.strand = +1; }
                    for( auto& iv : rev_intervals ) { iv.strand = -1; }
                    num_seeds_found = fwd_intervals.size() + rev_intervals.size();

                    // Pool both orientations into one peak_hits-sorted list before filtering,
                    // so SeedFilterConfig's thresholds (max_seeds, peak_hits_window, ...) apply
                    // globally across both strands rather than independently per strand.
                    spear::seeding::merge_seed_intervals(
                        fwd_intervals, rev_intervals, seed_intervals
                    );
                }

                // Discard candidates that don't meet the configured seed filter, before
                // spending alignment effort on them.
                spear::seeding::filter_seed_intervals( seed_intervals, seed_filter_cfg );
                num_seeds_kept = seed_intervals.size();

                // Empty result: read did not map (no k-mers, no valid interval, or every
                // candidate found got filtered out).
                if( seed_intervals.empty() ) {
                    if( keep_unmapped ) {
                        spear::align::AlignmentHit unmapped;
                        unmapped.ref_id  = -1;
                        unmapped.ref_pos = -1;
                        unmapped.flags   = spear::align::AlignmentFlags::Unmapped;
                        rr.hits.push_back( std::move( unmapped ));
                    }
                    bam_out_buff.emplace( read_idx, std::move( rr ));

                    if( emit_tsv ) {
                        line << file_base << '\t' << read_name << '\t' << read_length << '\t'
                             << num_kmers << '\t' << num_seeds_found << "\t0"
                             << "\t-\t-\t-\t-\t-\t-\t-\t-\t-\n";
                        tsv_out_buff->emplace( read_idx, line.str() );
                    }
                    return;
                }

                // Reverse complement of the read: needed for whichever candidates require the
                // opposite orientation (canonical mode always tries both against every candidate
                // window below; non-canonical mode uses it only for strand == -1 candidates).
                std::string const revcomp_seq = genesis::sequence::reverse_complement( seq );

                // Wfa2Aligner is not thread-safe, so each worker thread gets its own instance,
                // built once from the shared settings and reused across every read it processes.
                thread_local spear::align::Wfa2Aligner aligner( wfa2_settings );

                // Translate candidate intervals from genome-bin coordinates to sequence coordinates.
                for( auto& hit : seed_intervals ) {
                    auto const base_start = static_cast<std::uint64_t>( hit.left ) * w;
                    auto const base_end   = ( static_cast<std::uint64_t>( hit.right ) + 1 ) * w - 1;

                    auto const ref_seq = reference.find_by_global_position( base_start );
                    auto const seq_start = reference.global_offset_of( ref_seq.index );
                    auto const local_start = ref_seq.local_offset;
                    auto const local_end = std::min(
                        static_cast<std::size_t>( base_end ), seq_start + ref_seq.sequence->length - 1
                    ) - seq_start;

                    // Padded alignment window: the seed interval's own bounds only reflect where
                    // matching k-mers were *found* in the reference, not where in the read they
                    // came from, so the read's true alignment footprint can extend past
                    // [local_start, local_end] by up to a full read length on either side. Pad
                    // conservatively by read_length on each side, clamped so the window never
                    // goes negative or crosses into a neighboring contig.
                    std::size_t const win_start = ( local_start >= read_length )
                        ? ( local_start - read_length ) : 0;
                    std::size_t const win_end =
                        std::min( local_end + read_length, ref_seq.sequence->length - 1 );
                    // A view into the already-resident reference data, not a copy: Wfa2Aligner
                    // only needs the bytes for the duration of the call below.
                    std::string_view const target_window = std::string_view( ref_seq.sequence->sites )
                        .substr( win_start, win_end - win_start + 1 );

                    // Builds a BAM AlignmentHit from a successful Wfa2Result. ref_pos uses the
                    // alignment's own precise boundary (result.ref_begin, within target_window,
                    // translated back to a sequence-local coordinate via win_start) rather than
                    // the seed interval's coarser [local_start, local_end] bounds.
                    auto const make_alignment_hit_ = [&](
                        spear::align::Wfa2Result const& result, int strand
                    ) {
                        spear::align::AlignmentHit aln_hit;
                        aln_hit.ref_id  = static_cast<std::int32_t>( ref_seq.index );
                        aln_hit.ref_pos = static_cast<std::int64_t>( win_start ) + result.ref_begin;
                        aln_hit.flags   = ( strand < 0 ) ? spear::align::AlignmentFlags::Reverse : 0;
                        aln_hit.cigar   = result.cigar;
                        aln_hit.tag_as  = result.score;
                        aln_hit.tag_nm  = result.edit_distance;
                        return aln_hit;
                    };

                    // TSV diagnostics for this candidate; unused unless emit_tsv.
                    std::string align_score_str       = "-";
                    std::string align_score_other_str = "-";
                    std::string strand_tie_str        = "-";
                    std::string cigar_str             = "-";

                    if( canonical ) {
                        // Strand is not known from seeding alone here, so try both orientations
                        // against this candidate window and keep whichever scores better.
                        auto const fwd_result = aligner.align_cigar( seq, target_window );
                        auto const rev_result = aligner.align_cigar( revcomp_seq, target_window );
                        bool const fwd_ok =
                            fwd_result.status == spear::align::Wfa2Result::Status::Ok;
                        bool const rev_ok =
                            rev_result.status == spear::align::Wfa2Result::Status::Ok;

                        // Forward wins whenever only one side succeeded, and on exact ties.
                        bool const rev_wins =
                            rev_ok && ( !fwd_ok || rev_result.score > fwd_result.score );
                        bool const tie = fwd_ok && rev_ok && fwd_result.score == rev_result.score;

                        auto const& winner    = rev_wins ? rev_result : fwd_result;
                        auto const& loser     = rev_wins ? fwd_result : rev_result;
                        bool const winner_ok  = rev_wins ? rev_ok : fwd_ok;
                        bool const loser_ok   = rev_wins ? fwd_ok : rev_ok;

                        hit.strand = rev_wins ? -1 : +1;
                        if( winner_ok ) {
                            align_score_str = std::to_string( winner.score );
                            cigar_str       = spear::align::cigar_to_string( winner.cigar );
                            rr.hits.push_back( make_alignment_hit_( winner, hit.strand ));
                        }
                        if( loser_ok ) {
                            align_score_other_str = std::to_string( loser.score );
                        }
                        if( tie ) {
                            strand_tie_str = "tie";
                        }
                    } else {
                        // Strand is already known from which seeding query produced this
                        // candidate, so only one orientation needs to be tried here.
                        std::string_view const query_seq = ( hit.strand >= 0 )
                            ? std::string_view( seq ) : std::string_view( revcomp_seq );
                        auto const result = aligner.align_cigar( query_seq, target_window );
                        if( result.status == spear::align::Wfa2Result::Status::Ok ) {
                            align_score_str = std::to_string( result.score );
                            cigar_str       = spear::align::cigar_to_string( result.cigar );
                            rr.hits.push_back( make_alignment_hit_( result, hit.strand ));
                        }
                        // align_score_other_str and strand_tie_str stay "-": this path never
                        // tries the opposite orientation, so neither is applicable here.
                    }

                    if( emit_tsv ) {
                        // "0" should never appear here (strand is always set to +-1 above); printed
                        // literally, rather than defaulted to "+", as a visible check for that.
                        std::string const strand_str =
                            ( hit.strand > 0 ) ? "+" : ( hit.strand < 0 ) ? "-" : "0";

                        line << file_base << '\t' << read_name << '\t' << read_length << '\t'
                             << num_kmers << '\t' << num_seeds_found << '\t' << num_seeds_kept << '\t'
                             << ref_seq.sequence->name() << '\t'
                             << local_start << '\t' << local_end << '\t' << hit.peak_hits << '\t'
                             << strand_str << '\t'
                             << align_score_str << '\t' << align_score_other_str << '\t'
                             << strand_tie_str << '\t' << cigar_str << '\n'
                        ;
                    }
                }

                // Every candidate whose alignment succeeded is already in rr.hits; if none did
                // (all WFA2 attempts failed, e.g. score cap exceeded), this read ends up in the
                // same "no successful alignment" bucket as the empty-seed_intervals case above.
                if( rr.hits.empty() && keep_unmapped ) {
                    spear::align::AlignmentHit unmapped;
                    unmapped.ref_id  = -1;
                    unmapped.ref_pos = -1;
                    unmapped.flags   = spear::align::AlignmentFlags::Unmapped;
                    rr.hits.push_back( std::move( unmapped ));
                }
                bam_out_buff.emplace( read_idx, std::move( rr ));

                if( emit_tsv ) {
                    tsv_out_buff->emplace( read_idx, line.str() );
                }
            });
        }
    }

    // Wait for all reads to be fully processed before closing the output buffers.
    thread_pool->wait_for_all_pending_tasks();
    bam_out_buff.close();
    if( tsv_out_buff ) {
        tsv_out_buff->close();
    }

    // -------------------------------------------------------------------
    //     Statistics output for the user
    // -------------------------------------------------------------------

    LOG_MSG
        << "Processed " << total_reads << " read(s) from " << read_files.size()
        << " file(s)."
    ;
    LOG_MSG << seeding.stats().to_string();
    LOG_MSG << "Wrote results to " << options.output.get_output_filename( "map-align", bam_ext );
    if( options.emit_tsv.value ) {
        LOG_MSG << "Wrote diagnostic TSV to "
                << options.output.get_output_filename( "map-align", "tsv" );
    }
    global_options.report_memory_usage();
}

} // anonymous namespace

void run_map_align( MapAlignOptions const& options )
{
    using namespace genesis::util::format;
    using namespace genesis::util::io;

    // Read and parse the full JSON manifest.
    auto const& json_path = options.index_reader.index_file.value;
    auto const json_doc  = JsonReader().read( from_file( json_path ));
    auto const manifest  = manifest_from_json( json_doc );
    auto const sidx_path = resolve_sidx_path( json_path, manifest.metadata.sidx_path );

    // Reload the reference sequence data from the original FASTA files (the manifest only
    // stores names/lengths, not sites -- see spear::reference::build_collection()), using the
    // absolute, symlink-resolved paths recorded at index time. Validate the result against the
    // manifest's own metadata before using it: the genomic coordinates baked into the index
    // depend on the order and length of each sequence, so a stale or modified reference would
    // otherwise silently misalign every alignment against the wrong positions.

    LOG_MSG << "Loading reference genomes with " << manifest.reference.fasta_files().size() << " files.";
    auto reference = spear::reference::build_collection( manifest.reference.fasta_files(), true );
    if( !spear::reference::reference_collections_match( manifest.reference, reference ) ) {
        throw std::runtime_error(
            "Reference FASTA file(s) recorded in index '" + json_path + "' no longer match what "
            "the index was built from (different sequence count, order, name, or length). The "
            "reference may have moved, changed, or been replaced since `spear map index` was "
            "run; re-run `spear map index` against the current reference."
        );
    }
    LOG_MSG << "Loaded reference genomes with "
        << reference.size() << " sequence(s) totaling "
        << reference.total_length() << " bases."
    ;
    global_options.report_memory_usage();

    if( manifest.metadata.position_bits == 64 ) {
        run_map_align_impl<std::uint64_t>( options, manifest.metadata, reference, sidx_path );
    } else if( manifest.metadata.position_bits == 32 ) {
        run_map_align_impl<std::uint32_t>( options, manifest.metadata, reference, sidx_path );
    } else {
        throw std::runtime_error(
            "Index '" + json_path + "' has unsupported position_bits = "
            + std::to_string( manifest.metadata.position_bits )
        );
    }
}
