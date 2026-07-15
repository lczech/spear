#ifndef SPEAR_LIB_SEEDING_SEED_FILTER_H_
#define SPEAR_LIB_SEEDING_SEED_FILTER_H_

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
 * @brief Post-seeding filtering of candidate seed intervals, ahead of alignment.
 *
 * @file
 * @ingroup seeding
 */

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <vector>

namespace spear::seeding {

// =================================================================================================
//     SeedFilterConfig
// =================================================================================================

/**
 * @brief Configuration for filter_seed_intervals(): which candidate seed intervals from a
 * KmerSeeding query are worth pursuing further (e.g. aligning), out of everything found.
 *
 * All fields are independently optional; a disabled field imposes no constraint. Every seed
 * interval already satisfies peak_hits >= the min_hit_count used at seeding time (see
 * HitCollector::query()) -- these are additional, per-read-relative filters applied afterwards,
 * once the full set of candidates for a read (and hence its best peak_hits) is known.
 */
struct SeedFilterConfig
{
    /// Keep at most this many intervals (the strongest, since intervals are peak_hits-sorted
    /// descending). 0 = no limit.
    std::size_t max_seeds = 0;

    /// Discard intervals with peak_hits below this absolute value. 0 = disabled.
    std::size_t min_peak_hits = 0;

    /// Discard intervals whose peak_hits is more than this many below the best interval's
    /// peak_hits for the same read. Negative = disabled.
    int peak_hits_window = -1;

    /// Discard intervals whose peak_hits is below this fraction of the best interval's
    /// peak_hits for the same read. 0.0 = disabled.
    double peak_hits_fraction = 0.0;
};

// =================================================================================================
//     filter_seed_intervals
// =================================================================================================

/**
 * @brief Filter @p intervals in place, keeping only those that pass every active filter in
 * @p cfg, relative to the strongest interval in @p intervals.
 *
 * Requires @p intervals to already be sorted by peak_hits descending, as returned by
 * KmerSeeding::Query::finish_query() / HitCollector::query(). Every filter in SeedFilterConfig
 * is monotonic under that order: each one keeps some prefix of the list, so their combination
 * is also a prefix. This lets the implementation find the shortest surviving prefix with a
 * single forward scan and one resize(), rather than a general predicate-based removal.
 */
template<typename SeedIntervalT>
void filter_seed_intervals( std::vector<SeedIntervalT>& intervals, SeedFilterConfig const& cfg )
{
    if( intervals.empty() ) {
        return;
    }

    std::size_t const best = intervals.front().peak_hits;
    std::size_t cutoff = ( cfg.max_seeds > 0 )
        ? std::min( cfg.max_seeds, intervals.size() )
        : intervals.size();

    for( std::size_t i = 0; i < cutoff; ++i ) {
        auto const ph = intervals[i].peak_hits;
        bool const fails_min = cfg.min_peak_hits > 0 && ph < cfg.min_peak_hits;
        bool const fails_window = cfg.peak_hits_window >= 0 &&
            static_cast<long>( best - ph ) > cfg.peak_hits_window;
        bool const fails_fraction = cfg.peak_hits_fraction > 0.0 &&
            static_cast<double>( ph ) < cfg.peak_hits_fraction * static_cast<double>( best );
        if( fails_min || fails_window || fails_fraction ) {
            cutoff = i;
            break;
        }
    }

    intervals.resize( cutoff );
}

// =================================================================================================
//     merge_seed_intervals
// =================================================================================================

/**
 * @brief Merge two seed interval lists, each already sorted by peak_hits descending (e.g., one
 * per read orientation, from two separate KmerSeeding queries), into @p out, a single combined
 * list still sorted by peak_hits descending.
 *
 * Intended to run before filter_seed_intervals(), so that a single filtering pass applies its
 * thresholds (max_seeds, peak_hits_window, peak_hits_fraction, ...) globally across both inputs,
 * rather than independently per input.
 *
 * @p out is cleared first; it may alias neither @p a nor @p b.
 */
template<typename SeedIntervalT>
void merge_seed_intervals(
    std::vector<SeedIntervalT> const& a,
    std::vector<SeedIntervalT> const& b,
    std::vector<SeedIntervalT>&       out
) {
    auto const by_peak_hits_desc = []( SeedIntervalT const& x, SeedIntervalT const& y ) {
        return x.peak_hits > y.peak_hits;
    };
    assert(
        std::is_sorted( a.begin(), a.end(), by_peak_hits_desc )
        && "merge_seed_intervals: 'a' must be sorted by peak_hits descending"
    );
    assert(
        std::is_sorted( b.begin(), b.end(), by_peak_hits_desc )
        && "merge_seed_intervals: 'b' must be sorted by peak_hits descending"
    );

    out.clear();
    out.reserve( a.size() + b.size() );
    std::merge( a.begin(), a.end(), b.begin(), b.end(), std::back_inserter( out ), by_peak_hits_desc );
}

/**
 * @brief Convenience overload returning results by value.
 *
 * Prefer the out-parameter overload with a reused buffer to avoid per-call allocation on the
 * hot path.
 */
template<typename SeedIntervalT>
std::vector<SeedIntervalT> merge_seed_intervals(
    std::vector<SeedIntervalT> const& a,
    std::vector<SeedIntervalT> const& b
) {
    std::vector<SeedIntervalT> out;
    merge_seed_intervals( a, b, out );
    return out;
}

} // namespace spear::seeding

#endif // SPEAR_LIB_SEEDING_SEED_FILTER_H_
