#ifndef SPEAR_LIB_SEEDING_KMER_SEEDING_H_
#define SPEAR_LIB_SEEDING_KMER_SEEDING_H_

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

/**
 * @brief K-mer seeding engine over an InvertedIndex.
 *
 * @file
 * @ingroup seeding
 */

#include "spear/inverted_index/hit_collector.hpp"
#include "spear/inverted_index/index.hpp"
#include "spear/inverted_index/term_postings.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace spear::seeding {

// =================================================================================================
//     KmerSeeding
// =================================================================================================

/**
 * @brief Finds candidate genomic intervals for a read by k-mer seeding against an InvertedIndex.
 *
 * Abstracts TermPostings and HitCollector behind genomics/k-mer terminology. Thread-safe:
 * multiple threads may call start_query() concurrently; per-thread state (TermPostings,
 * HitCollector) is managed as thread-local storage internally. Each call to start_query()
 * returns a Query object that drives one seeding round via add_kmer() for a read, then call
 * finish_query().
 *
 * Construct once per index, then call start_query() from as many threads as needed:
 * @code
 *   KmerSeeding<uint64_t> seeding( config, index, bin_width );
 *   // In each thread:
 *   auto q = seeding.start_query();
 *   for( auto const& kmer : extractor ) {
 *       if( !q.add_kmer( term_index ) ) break;   // optional early exit on truncation
 *   }
 *   std::vector<KmerSeeding<uint64_t>::SeedInterval> out;
 *   q.finish_query( read_length, out );
 * @endcode
 *
 * @tparam PositionT  Unsigned integer type matching the InvertedIndex; uint32_t or uint64_t.
 */
template<typename PositionT = std::uint64_t>
class KmerSeeding
{
    static_assert(
        std::is_same_v<PositionT, std::uint32_t> || std::is_same_v<PositionT, std::uint64_t>,
        "KmerSeeding requires PositionT = uint32_t or uint64_t"
    );

public:

    // -------------------------------------------------------------------------
    //     Constants
    // -------------------------------------------------------------------------

    // Ring buffer capacity for the internal HitCollector.
    // Supports window lengths up to kRingCap - 1 genome bins.
    static constexpr std::size_t kRingCap = 32;

    // Maximum k-mers per read. Reads with more k-mers are seeded using only the first kMaxKmersPerRead;
    // the truncation is recorded in stats().truncated_queries for the caller to report.
    static constexpr std::size_t kMaxKmersPerRead =
        spear::inverted_index::HitCollector<PositionT, kRingCap>::kMaxLists;

    // Maximum window length (in genome-bin units) the internal HitCollector can handle.
    static constexpr std::size_t kMaxWindowLength = kRingCap - 1;

    // -------------------------------------------------------------------------
    //     Types
    // -------------------------------------------------------------------------

    // Result interval type from one seeding query. Fields are in genome-bin coordinates;
    // translate to reference-sequence base coordinates at the call site as needed.
    using SeedInterval = typename spear::inverted_index::HitCollector<PositionT, kRingCap>::HitInterval;

    // Encoded k-mer index type, matching InvertedIndex::term_index_type.
    using kmer_index_type = typename spear::inverted_index::InvertedIndex<PositionT>::term_index_type;

    // -------------------------------------------------------------------------
    //     Config
    // -------------------------------------------------------------------------

    /**
     * @brief Configuration for KmerSeeding, typically populated from CLI options.
     *
     * All fields have defaults that produce reasonable seeding behaviour without tuning:
     * - min_hit_count and min_hit_fraction both 0 → require at least 1 k-mer hit per window.
     * - max_occurrences_per_kmer 0 → no soft cap; all k-mers are looked up.
     */
    struct Config
    {
        /// Minimum number of distinct k-mer hits required within a window.
        /// 0 means use min_hit_fraction if set, otherwise default to 1.
        std::size_t min_hit_count            = 0;

        /// Fraction of a read's k-mers required within a window. Used when min_hit_count == 0.
        /// 0.0 means use the default of 1.
        double      min_hit_fraction         = 0.0;

        /// Soft cap on posting-list length: k-mers occurring more times than this are skipped.
        /// 0 means no cap.
        std::size_t max_occurrences_per_kmer = 0;

        /// Whether finish_query() accumulates stats() bookkeeping. Disable only to measure the
        /// overhead of stats collection itself, or to shave it off on a performance-critical run.
        bool collect_stats = true;
    };

    // -------------------------------------------------------------------------
    //     Stats
    // -------------------------------------------------------------------------

    /**
     * @brief Cumulative statistics across all queries run through this engine.
     *
     * Plain value snapshot returned by stats(). Safe to copy, print, and pass around freely.
     */
    struct Stats
    {
        /// Total number of finish_query() calls.
        std::uint64_t total_queries     = 0;

        /// Queries where the read had more than kMaxKmersPerRead k-mers;
        /// only the first kMaxKmersPerRead were seeded.
        std::uint64_t truncated_queries = 0;

        /// Queries that returned no seed intervals (read did not map).
        std::uint64_t empty_queries     = 0;

        /// Cumulative k-mers whose posting list was found (present, decoded) across all queries.
        std::uint64_t found_kmers       = 0;

        /// Cumulative k-mers not present in the index at all (zero postings) across all queries.
        std::uint64_t absent_kmers      = 0;

        /// Cumulative k-mers whose posting list exceeded the index's build-time hard cap
        /// (a property of the index itself, unlike soft_capped_kmers which is tunable per run).
        std::uint64_t hard_capped_kmers = 0;

        /// Cumulative k-mers skipped due to max_occurrences_per_kmer across all queries.
        std::uint64_t soft_capped_kmers = 0;

        /// Count of reads by number of k-mers added to the query (bin i = reads with exactly
        /// i k-mers added), for i in [0, kMaxKmersPerRead]. Always has kMaxKmersPerRead + 1
        /// entries, since the bound is a compile-time constant.
        std::vector<std::uint64_t> num_kmers_histogram;

        /// Count of reads by the peak_hits of their single strongest seed interval (bin i =
        /// reads whose best interval had peak_hits == i), for i in [1, kMaxKmersPerRead]. Bin 0
        /// counts reads with no seed interval at all: a real interval always has peak_hits >= 1
        /// (HitCollector::query() never returns one below min_hit_count), so bin 0 is otherwise
        /// unreachable and is repurposed here as the "read did not map" bucket.
        std::vector<std::uint64_t> top_peak_hits_histogram;

        /// Count of reads by number of seed intervals returned (bin i = reads with exactly i
        /// candidate intervals). Unbounded: a read can map to arbitrarily many intervals in
        /// repetitive genome regions, so this grows as large as the largest count observed
        /// across all queries so far.
        std::vector<std::uint64_t> num_seeds_histogram;

        /**
         * @brief Human-readable dump of every field of the stats.
         *
         * Scalar counters are always printed, including zeros. Histogram bins that are zero
         * are skipped (a display compaction, not a hidden field: every histogram is still
         * present, just without empty entries).
         */
        std::string to_string() const
        {
            std::ostringstream os;
            os << "Total queries: " << total_queries << '\n';
            os << "Truncated queries (more than " << KmerSeeding::kMaxKmersPerRead
               << " k-mers): " << truncated_queries << '\n';
            os << "Empty queries (no seed interval found): " << empty_queries << '\n';
            os << "K-mers found in the index: " << found_kmers << '\n';
            os << "K-mers not present in the index at all: " << absent_kmers << '\n';
            os << "K-mers skipped due to hard cap: " << hard_capped_kmers << '\n';
            os << "K-mers skipped due to soft cap: " << soft_capped_kmers << '\n';

            auto const print_histogram = [&](
                char const* label, std::vector<std::uint64_t> const& hist
            ) {
                os << label << '\n';
                for( std::size_t i = 0; i < hist.size(); ++i ) {
                    if( hist[i] == 0 ) {
                        continue;
                    }
                    os << "  " << i << ": " << hist[i] << '\n';
                }
            };
            print_histogram( "Histogram of k-mers per read:", num_kmers_histogram );
            print_histogram( "Histogram of candidate intervals per read:", num_seeds_histogram );
            print_histogram(
                "Histogram of top seed interval strength (peak_hits) per read "
                "(0 = read did not map):",
                top_peak_hits_histogram
            );

            return os.str();
        }
    };

    // -------------------------------------------------------------------------
    //     Query
    // -------------------------------------------------------------------------

    /**
     * @brief Per-read seeding handle returned by start_query().
     *
     * Call add_kmer() for each encoded k-mer in the read, then finish_query() to get the
     * seed intervals. The Query must not be destroyed before finish_query() is called —
     * doing so calls std::terminate() to catch the programming error loudly.
     *
     * Move-only. The moved-from Query is considered finished and its destructor is a no-op.
     */
    class Query
    {
    public:
        Query( Query const& ) = delete;
        Query& operator=( Query const& ) = delete;
        Query& operator=( Query&& ) = delete;

        Query( Query&& other ) noexcept
            : engine_( other.engine_ )
            , term_postings_( other.term_postings_ )
            , hit_collector_( other.hit_collector_ )
            , index_( other.index_ )
            , count_( other.count_ )
            , truncated_( other.truncated_ )
            , finished_( other.finished_ )
        {
            // Moved-from object must not call terminate in its destructor
            other.finished_ = true;
        }

        ~Query() noexcept
        {
            if( !finished_ ) {
                std::fprintf( stderr,
                    "KmerSeeding::Query destroyed without calling finish_query() — "
                    "per-thread seeding state is now corrupt\n"
                );
                std::terminate();
            }
        }

        /**
         * @brief Add one k-mer (by encoded term index) to this seeding round.
         *
         * @returns true if the k-mer was accepted; false if the read already has
         * kMaxKmersPerRead k-mers and the k-mer was silently dropped. On the first
         * dropped call the truncation is flagged and reported in stats().truncated_queries.
         * The caller may break the extraction loop early on false if desired.
         */
        bool add_kmer( KmerSeeding::kmer_index_type term_index )
        {
            if( count_ >= KmerSeeding::kMaxKmersPerRead ) {
                truncated_ = true;
                return false;
            }
            term_postings_->add_deferred( index_, term_index );
            ++count_;
            return true;
        }

        /**
         * @brief Flush, run the hit-collector query, write results to @p out, and commit stats.
         *
         * @param read_length  Length of the read in bases; used to auto-derive the window length.
         * @param out          Output vector; cleared then filled with seed intervals sorted by
         *                     peak_hits descending.
         *
         * Must be called exactly once per Query. Calling it a second time, or not calling it
         * before the Query is destroyed, is a programming error.
         */
        void finish_query( std::size_t read_length, std::vector<SeedInterval>& out )
        {
            assert( !finished_ );

            // Mark finished and clear the per-thread in-use guard before any throwing code,
            // so that start_query() can be called again even if an exception propagates out.
            finished_ = true;
            in_use_flag_() = false;

            // Flush deferred posting-list lookups (prefetch-pipelined)
            term_postings_->flush_deferred( index_ );

            // Derive per-read parameters and run the hit collector
            auto const window_length = engine_.derive_window_length_( read_length );
            auto const min_hits      = engine_.derive_min_hits_( count_ );
            hit_collector_->query( *term_postings_, window_length, min_hits, out );

            // Stats collection can be switched off (Config::collect_stats); seeding itself
            // (out, filled above) is unaffected either way.
            if( !engine_.config_.collect_stats ) {
                return;
            }

            // Accumulate stats
            engine_.stats_.total_queries.fetch_add( 1, std::memory_order_relaxed );
            if( truncated_ ) {
                engine_.stats_.truncated_queries.fetch_add( 1, std::memory_order_relaxed );
            }
            if( out.empty() ) {
                engine_.stats_.empty_queries.fetch_add( 1, std::memory_order_relaxed );
            }
            auto const& tp_stats = term_postings_->stats();
            engine_.stats_.hard_capped_kmers.fetch_add(
                static_cast<std::uint64_t>( tp_stats.hard_capped ), std::memory_order_relaxed
            );
            engine_.stats_.soft_capped_kmers.fetch_add(
                static_cast<std::uint64_t>( tp_stats.soft_capped ), std::memory_order_relaxed
            );
            engine_.stats_.found_kmers.fetch_add(
                static_cast<std::uint64_t>( tp_stats.found ), std::memory_order_relaxed
            );
            engine_.stats_.absent_kmers.fetch_add(
                static_cast<std::uint64_t>( tp_stats.empty ), std::memory_order_relaxed
            );

            // Bounded histograms: count_ <= kMaxKmersPerRead is guaranteed by add_kmer()'s own
            // limit, and peak_hits of any returned interval is a popcount over a kMaxKmersPerRead-
            // bit set, so both indices are in range by construction; the asserts are a cheap
            // safety net against that invariant ever breaking under future changes.
            assert( count_ < engine_.stats_.num_kmers_histogram.size() );
            engine_.stats_.num_kmers_histogram[count_].fetch_add( 1, std::memory_order_relaxed );

            // peak_hits of the strongest (first, since out is sorted descending) interval; bin 0
            // is reserved for "no seed interval found" (see Stats::top_peak_hits_histogram).
            std::size_t const top_peak_hits = out.empty() ? 0 : out[0].peak_hits;
            assert( top_peak_hits < engine_.stats_.top_peak_hits_histogram.size() );
            engine_.stats_.top_peak_hits_histogram[top_peak_hits].fetch_add(
                1, std::memory_order_relaxed
            );

            // Unbounded histogram of seed-interval counts: growth can't be done lock-free, so
            // this is the one stat guarded by a mutex instead of an atomic. One lock per read,
            // negligible next to the k-mer lookups and sliding-window merge just performed above.
            {
                std::lock_guard<std::mutex> lock( engine_.stats_.num_seeds_histogram_mutex );
                auto& hist = engine_.stats_.num_seeds_histogram;
                if( out.size() >= hist.size() ) {
                    hist.resize( out.size() + 1, 0 );
                }
                ++hist[ out.size() ];
            }
        }

        /**
         * @brief Convenience overload returning results by value.
         *
         * Prefer the out-parameter overload with a reused thread_local vector to avoid
         * per-query allocation on the hot path.
         */
        std::vector<SeedInterval> finish_query( std::size_t read_length )
        {
            std::vector<SeedInterval> out;
            finish_query( read_length, out );
            return out;
        }

    private:

        friend class KmerSeeding;

        Query(
            KmerSeeding const&                                          engine,
            spear::inverted_index::TermPostings<PositionT>&             term_postings,
            spear::inverted_index::HitCollector<PositionT, kRingCap>&   hit_collector,
            spear::inverted_index::InvertedIndex<PositionT> const&      index
        )
            : engine_( engine )
            , term_postings_( &term_postings )
            , hit_collector_( &hit_collector )
            , index_( index )
        {}

        // Per-thread in-use guard: set by start_query(), cleared by finish_query().
        // Separate per PositionT instantiation so uint32_t and uint64_t engines
        // each have their own flag and do not conflict.
        static bool& in_use_flag_() noexcept
        {
            thread_local bool flag = false;
            return flag;
        }

        KmerSeeding const&                                              engine_;
        spear::inverted_index::TermPostings<PositionT>*                 term_postings_;
        spear::inverted_index::HitCollector<PositionT, kRingCap>*       hit_collector_;
        spear::inverted_index::InvertedIndex<PositionT> const&          index_;
        std::size_t count_     = 0;
        bool        truncated_ = false;
        bool        finished_  = false;
    };

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    /**
     * @brief Construct a KmerSeeding engine.
     *
     * @param config     Tuning parameters (min hits, soft cap).
     * @param index      Opened InvertedIndex; must outlive this engine.
     * @param bin_width  Genome bin width (in bases) used when the index was built; used to
     *                   auto-derive the window length per read.
     *
     * @throws std::invalid_argument if bin_width == 0.
     */
    KmerSeeding(
        Config const&                                               config,
        spear::inverted_index::InvertedIndex<PositionT> const&      index,
        std::size_t                                                 bin_width
    )
        : config_( config )
        , index_( index )
        , bin_width_( bin_width )
    {
        if( bin_width_ == 0 ) {
            throw std::invalid_argument( "KmerSeeding: bin_width must be > 0" );
        }
    }

    // Non-copyable and non-moveable: holds a reference to the index and owns non-moveable Stats.
    KmerSeeding( KmerSeeding const& ) = delete;
    KmerSeeding( KmerSeeding&& )      = delete;
    KmerSeeding& operator=( KmerSeeding const& ) = delete;
    KmerSeeding& operator=( KmerSeeding&& )      = delete;

    // -------------------------------------------------------------------------
    //     Query interface
    // -------------------------------------------------------------------------

    /**
     * @brief Start a new seeding round for one read.
     *
     * Prepares the thread-local TermPostings and HitCollector for a new query and returns
     * a Query object. Only one Query may be active per thread at a time; starting a second
     * one before calling finish_query() on the first throws std::runtime_error.
     *
     * @returns A [[nodiscard]] Query; must be captured and finish_query() must be called.
     */
    [[nodiscard]] Query start_query() const
    {
        if( Query::in_use_flag_() ) {
            throw std::runtime_error(
                "KmerSeeding::start_query: a Query is already active on this thread; "
                "call finish_query() before starting a new one"
            );
        }
        Query::in_use_flag_() = true;

        // Thread-local TermPostings and HitCollector: allocated once per thread on the
        // first call, then reused across all subsequent queries on the same thread.
        thread_local spear::inverted_index::TermPostings<PositionT> term_postings( kMaxKmersPerRead );
        thread_local spear::inverted_index::HitCollector<PositionT, kRingCap> hit_collector;

        // Reset for the new round (keeps buffer allocations intact)
        term_postings.clear();
        term_postings.set_max_posting_list_length(
            static_cast<std::uint64_t>( config_.max_occurrences_per_kmer )
        );

        return Query( *this, term_postings, hit_collector, index_ );
    }

    /**
     * @brief Run one seeding query in a single call, for callers that already hold all of a
     * read's k-mers in a contiguous container (as opposed to streaming them one at a time).
     *
     * Equivalent to calling start_query(), add_kmer() for each element of @p kmers, then
     * finish_query( read_length, out ). If @p kmers exceeds kMaxKmersPerRead, seeding uses
     * only the first kMaxKmersPerRead entries and the truncation is recorded in stats(), same
     * as the streaming API.
     *
     * @param kmers        Encoded k-mer indices for the read, in order.
     * @param read_length  Length of the read in bases; used to auto-derive the window length.
     * @param out          Output vector; cleared then filled with seed intervals sorted by
     *                     peak_hits descending.
     */
    void run_query(
        std::span<kmer_index_type const> kmers,
        std::size_t                      read_length,
        std::vector<SeedInterval>&       out
    ) const {
        auto q = start_query();
        for( auto const kmer : kmers ) {
            if( !q.add_kmer( kmer ) ) {
                break; // truncated_ already flagged by add_kmer; further calls are no-ops
            }
        }
        q.finish_query( read_length, out );
    }

    /**
     * @brief Convenience overload of run_query() returning results by value.
     *
     * Prefer the out-parameter overload with a reused thread_local vector to avoid
     * per-query allocation on the hot path.
     */
    std::vector<SeedInterval> run_query( std::span<kmer_index_type const> kmers, std::size_t read_length ) const
    {
        std::vector<SeedInterval> out;
        run_query( kmers, read_length, out );
        return out;
    }

    // -------------------------------------------------------------------------
    //     Stats
    // -------------------------------------------------------------------------

    /**
     * @brief Snapshot of cumulative statistics across all finish_query() calls on this engine.
     *
     * Only call this once every finish_query() call that was started has also returned, on
     * every thread that used this engine -- e.g. after joining/waiting on whatever thread pool
     * drove the queries. Calling it while a query is still in flight elsewhere is not safe.
     */
    Stats stats() const
    {
        Stats s;
        s.total_queries     = stats_.total_queries.load( std::memory_order_relaxed );
        s.truncated_queries = stats_.truncated_queries.load( std::memory_order_relaxed );
        s.empty_queries     = stats_.empty_queries.load( std::memory_order_relaxed );
        s.soft_capped_kmers = stats_.soft_capped_kmers.load( std::memory_order_relaxed );
        s.found_kmers       = stats_.found_kmers.load( std::memory_order_relaxed );
        s.absent_kmers      = stats_.absent_kmers.load( std::memory_order_relaxed );
        s.hard_capped_kmers = stats_.hard_capped_kmers.load( std::memory_order_relaxed );

        s.num_kmers_histogram.reserve( stats_.num_kmers_histogram.size() );
        for( auto const& bin : stats_.num_kmers_histogram ) {
            s.num_kmers_histogram.push_back( bin.load( std::memory_order_relaxed ) );
        }

        s.top_peak_hits_histogram.reserve( stats_.top_peak_hits_histogram.size() );
        for( auto const& bin : stats_.top_peak_hits_histogram ) {
            s.top_peak_hits_histogram.push_back( bin.load( std::memory_order_relaxed ) );
        }

        {
            std::lock_guard<std::mutex> lock( stats_.num_seeds_histogram_mutex );
            s.num_seeds_histogram = stats_.num_seeds_histogram;
        }

        return s;
    }

    // -------------------------------------------------------------------------
    //     Private helpers
    // -------------------------------------------------------------------------

private:

    // Derive window length in genome-bin units from read_length (bases) and bin_width_.
    // Auto-derives as ceil(read_length / bin_width) + 1, clamped to [1, kMaxWindowLength].
    PositionT derive_window_length_( std::size_t read_length ) const
    {
        std::size_t const raw = static_cast<std::size_t>(
            std::ceil( static_cast<double>( read_length ) / static_cast<double>( bin_width_ ))
        ) + 1;
        return static_cast<PositionT>( std::clamp( raw, std::size_t{ 1 }, kMaxWindowLength ) );
    }

    // Derive minimum hit count from the number of k-mers added in this query round.
    // Priority: explicit min_hit_count > fraction-derived > default of 1.
    std::size_t derive_min_hits_( std::size_t num_kmers ) const
    {
        if( config_.min_hit_count > 0 ) {
            return config_.min_hit_count;
        }
        if( config_.min_hit_fraction > 0.0 ) {
            auto const m = static_cast<std::size_t>(
                std::llround( static_cast<double>( num_kmers ) * config_.min_hit_fraction )
            );
            return std::max( m, std::size_t{ 1 } );
        }
        return 1;
    }

    // -------------------------------------------------------------------------
    //     Data members
    // -------------------------------------------------------------------------

    Config config_;
    spear::inverted_index::InvertedIndex<PositionT> const& index_;
    std::size_t bin_width_;

    // Mutable state backing Stats, written from finish_query() and read (snapshotted) by
    // stats(). Kept internal so callers never have to deal with atomics or the histogram lock
    // directly (can't copy/print/return an atomic by value, and the plain scalars/histograms
    // in Stats need somewhere thread-safe to live). Most fields here are lock-free atomics;
    // num_seeds_histogram is the one exception (see its own comment below).
    struct StatsStorage
    {
        std::atomic<std::uint64_t> total_queries     = 0;
        std::atomic<std::uint64_t> truncated_queries = 0;
        std::atomic<std::uint64_t> empty_queries     = 0;
        std::atomic<std::uint64_t> found_kmers       = 0;
        std::atomic<std::uint64_t> absent_kmers      = 0;
        std::atomic<std::uint64_t> hard_capped_kmers = 0;
        std::atomic<std::uint64_t> soft_capped_kmers = 0;

        // Bounded histograms: index is always in range by construction (see finish_query()),
        // so a plain fixed-size array with per-element atomic increments is enough -- no lock
        // needed, unlike num_seeds_histogram below.
        std::array<std::atomic<std::uint64_t>, kMaxKmersPerRead + 1> num_kmers_histogram{};
        std::array<std::atomic<std::uint64_t>, kMaxKmersPerRead + 1> top_peak_hits_histogram{};

        // Histogram of seed-interval counts (Stats::num_seeds_histogram). Unlike the bounded
        // histograms above, the number of seed intervals for a read has no fixed upper bound,
        // so this can't be a fixed-size array of atomics (std::atomic isn't movable, so a
        // std::vector of them can't grow); growing a plain vector isn't lock-free either, so
        // this is instead a plain vector guarded by its own mutex (see finish_query() and
        // stats()), the one field here that isn't a lock-free atomic.
        std::mutex num_seeds_histogram_mutex;
        std::vector<std::uint64_t> num_seeds_histogram;
    };

    // Mutable so Query (which holds KmerSeeding const&) can update this in finish_query().
    mutable StatsStorage stats_;
};

} // namespace spear::seeding

#endif // SPEAR_LIB_SEEDING_KMER_SEEDING_H_
