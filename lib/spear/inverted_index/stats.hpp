#ifndef SPEAR_LIB_INVERTED_INDEX_STATS_H_
#define SPEAR_LIB_INVERTED_INDEX_STATS_H_

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
 * @brief Statistics structs and accumulation helpers for the Inverted Index components.
 *
 * @file
 * @ingroup util
 */

#include <cstddef>
#include <cstdint>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Inverted Index Stats
// =================================================================================================

/**
 * @brief Statistics produced by InvertedIndexBuilder::write().
 *
 * The posting_count_histogram has size max_postings_per_term + 2:
 *   - histogram[0]        : number of terms with no postings (empty)
 *   - histogram[1..max]   : number of terms with exactly that many postings
 *   - histogram[max + 1]  : number of terms that were capped
 *
 * total_blob_bytes is the total compressed size of all posting lists,
 * not counting the offset table or footer.
 */
struct InvertedIndexStats
{
    std::vector<std::uint64_t> posting_count_histogram;
    std::uint64_t              total_blob_bytes = 0;
};

// =================================================================================================
//     Hit Collector Stats
// =================================================================================================

/**
 * @brief Accumulated statistics across one or more HitCollector::query() calls.
 *
 * All histograms grow on demand to fit the largest observed value; iterate up to `.size()`
 * to print or process them without needing to know the maximum in advance.
 *
 * Per-interval histograms (one entry per interval across all queries):
 *
 *   interval_length_histogram[n]  : number of intervals with right - left == n
 *   peak_count_histogram[n]       : number of intervals whose peak distinct-list count == n
 *
 * Per-query histograms (one entry per accumulate() call, e.g. per read):
 *
 *   max_peak_histogram[n]              : number of queries whose best interval had peak_hits == n;
 *                                        n == 0 means the query returned no intervals at all
 *                                        (the read produced no hit above the min_hit_count threshold)
 *   intervals_per_query_histogram[n]   : number of queries that returned exactly n intervals;
 *                                        n == 0 counts unmapped reads directly
 *   total_evidence_histogram[n]        : number of queries whose sum of peak_hits across all
 *                                        intervals equals n; distinguishes a read with one strong
 *                                        hit (low sum, high max_peak) from one with many weak hits
 *                                        (high sum, low max_peak) typical of repetitive regions
 *
 * Scalar counters:
 *
 *   total_intervals  : total intervals emitted across all queries
 *   total_queries    : number of times accumulate() was called (i.e. number of reads queried)
 */
struct HitCollectorStats
{
    // Per-interval
    std::vector<std::uint64_t> interval_length_histogram;
    std::vector<std::uint64_t> peak_count_histogram;

    // Per-query (per-read)
    std::vector<std::uint64_t> max_peak_histogram;
    std::vector<std::uint64_t> intervals_per_query_histogram;
    std::vector<std::uint64_t> total_evidence_histogram;

    // Scalar totals
    std::uint64_t total_intervals = 0;
    std::uint64_t total_queries   = 0;
};

// =================================================================================================
//     Hit Collector Stats Accumulation
// =================================================================================================

/**
 * @brief Add the results of one HitCollector::query() call to @p stats.
 *
 * Increments total_queries by one regardless of whether any intervals were found.
 * Updates all four histograms and total_intervals.
 *
 * Per-interval: interval_length_histogram and peak_count_histogram are updated for every interval.
 * Per-query: intervals_per_query_histogram is updated with the number of intervals returned;
 * max_peak_histogram is updated with the highest peak_hits seen (0 if no intervals).
 *
 * @p intervals is typically the output vector from HitCollector::query(); any type whose
 * elements expose `left`, `right`, and `peak_hits` fields is accepted.
 */
template<typename HitInterval>
void accumulate( HitCollectorStats& stats, std::vector<HitInterval> const& intervals )
{
    // Helper: resize histogram on demand and increment the bucket at index.
    auto tally = []( std::vector<std::uint64_t>& hist, std::size_t index )
    {
        if( index >= hist.size() ) {
            hist.resize( index + 1, 0 );
        }
        ++hist[index];
    };

    // Totals.
    ++stats.total_queries;

    // Per-interval pass; track per-query aggregates on the fly to avoid a second scan.
    std::size_t max_peak       = 0;
    std::size_t total_evidence = 0;
    for( auto const& iv : intervals ) {
        auto const pc = static_cast<std::size_t>( iv.peak_hits );
        if( pc > max_peak ) {
            max_peak = pc;
        }
        total_evidence += pc;
        tally( stats.interval_length_histogram, static_cast<std::size_t>( iv.right - iv.left ));
        tally( stats.peak_count_histogram, pc );
        ++stats.total_intervals;
    }

    // Per-query summaries (all stay 0 when intervals is empty, which is correct)
    tally( stats.max_peak_histogram,            max_peak          );
    tally( stats.intervals_per_query_histogram, intervals.size()  );
    tally( stats.total_evidence_histogram,      total_evidence    );
}

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_STATS_H_
