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
 * Histograms grow on demand to fit the largest observed value; iterate up to `.size()`
 * to print or process them without needing to know the maximum in advance.
 *
 *   interval_length_histogram[n]  : number of intervals with v_right - v_left == n
 *   peak_count_histogram[n]       : number of intervals whose peak distinct-list count == n
 *
 * total_intervals counts every interval emitted across all accumulated queries.
 * total_queries counts the number of times accumulate() has been called, allowing
 * per-query averages to be derived.
 */
struct HitCollectorStats
{
    std::vector<std::uint64_t> interval_length_histogram;
    std::vector<std::uint64_t> peak_count_histogram;
    std::uint64_t              total_intervals = 0;
    std::uint64_t              total_queries   = 0;
};

// =================================================================================================
//     Hit Collector Stats Accumulation
// =================================================================================================

/**
 * @brief Add the results of one HitCollector::query() call to @p stats.
 *
 * Increments total_queries by one regardless of whether any intervals were found.
 * For each interval, increments total_intervals and extends and updates both histograms.
 *
 * @p intervals is typically the output vector from HitCollector::query(); any type whose
 * elements expose `v_left`, `v_right`, and `peak_count` fields is accepted.
 */
template<typename HitInterval>
void accumulate( HitCollectorStats& stats, std::vector<HitInterval> const& intervals )
{
    ++stats.total_queries;
    for( auto const& iv : intervals ) {
        ++stats.total_intervals;

        auto const len = static_cast<std::size_t>( iv.v_right - iv.v_left );
        if( len >= stats.interval_length_histogram.size() ) {
            stats.interval_length_histogram.resize( len + 1, 0 );
        }
        ++stats.interval_length_histogram[len];

        auto const pc = static_cast<std::size_t>( iv.peak_count );
        if( pc >= stats.peak_count_histogram.size() ) {
            stats.peak_count_histogram.resize( pc + 1, 0 );
        }
        ++stats.peak_count_histogram[pc];
    }
}

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_STATS_H_
