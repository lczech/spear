#ifndef SPEAR_LIB_INVERTED_INDEX_HIT_COLLECTOR_H_
#define SPEAR_LIB_INVERTED_INDEX_HIT_COLLECTOR_H_

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
 * @brief Streaming multi-term hit collection and interval detection over an InvertedIndex.
 *
 * @file
 * @ingroup util
 */

#include "spear/inverted_index/index.hpp"
#include "spear/inverted_index/term_postings.hpp"
#include "genesis/util/container/tournament_tree.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Hit Collector
// =================================================================================================

/**
 * @brief Finds maximal hit intervals across a TermPostings collection.
 *
 * Given a TermPostings collection, a window length, and a minimum hit count,
 * finds every maximal interval [v_left, v_right] in posting-value space such that a
 * sliding window of size window_length can be placed anywhere within it and still see at least
 * min_hit_count distinct posting lists represented.
 *
 * Internally runs a streaming k-way merge over the S sorted posting lists via a
 * genesis::util::container::TournamentTree, maintaining a small ring buffer (size @p RingCap,
 * a power of two, such that RingCap > window_length) of (value, ListBitset) window entries.
 * All scratch state (tournament tree, ring buffer) is held as members and reused across calls
 * to avoid per-query allocations.
 *
 * @tparam PositionT  Unsigned integer type for posting positions. numeric_limits<PositionT>::max()
 *                    is reserved as a sentinel by the tournament tree and must not appear as an
 *                    actual posting value; InvertedIndex::open() enforces this via max_position().
 * @tparam RingCap    Ring buffer capacity; must be a power of two and satisfy
 *                    RingCap > window_length for every window_length passed to query().
 *                    Default 32 supports window_length up to 31; use 64 for up to 63, etc.
 */
template<typename PositionT = std::uint64_t, std::size_t RingCap = 32>
class HitCollector
{
    static_assert( (RingCap & (RingCap - 1)) == 0, "RingCap must be a power of two" );

public:

    // -------------------------------------------------------------------------
    //     Constants
    // -------------------------------------------------------------------------

    // Maximum number of distinct posting lists (k-mers) trackable per query window.
    // ListBitset derives its word-array size from this; changing it here propagates automatically.
    static constexpr std::size_t kMaxLists = 256;
    static_assert( kMaxLists % 64 == 0, "kMaxLists must be a multiple of 64" );

    // -------------------------------------------------------------------------
    //     Typedefs
    // -------------------------------------------------------------------------

    using position_type = PositionT;

    /**
     * @brief A maximal interval [v_left, v_right] where the hit condition is continuously met.
     *
     * Every window of size window_length placed within this interval captures at least
     * min_hit_count distinct lists.
     * `peak_count` is the maximum distinct-list count seen across all window positions within
     * the interval; used to sort results by strength (highest peak first).
     */
    struct HitInterval
    {
        PositionT   v_left;
        PositionT   v_right;
        std::size_t peak_count = 0;
    };

    // -------------------------------------------------------------------------
    //     Query
    // -------------------------------------------------------------------------

    /**
     * @brief Find all maximal hit intervals; results are appended to @p hit_intervals
     * (cleared first), and sorted by peak_count descending, then v_left ascending.
     *
     * @param term_postings   Populated TermPostings collection.
     * @param window_length   Window length: a window [x, x+window_length] must cover values
     *                        from >= min_hit_count lists.
     * @param min_hit_count   Minimum number of distinct lists required; must be >= 1.
     * @param hit_intervals   Output vector; cleared then filled with results. Sorted by peak
     *                        hit count.
     *
     * Returns immediately with an empty @p hit_intervals if min_hit_count > list_count().
     * Throws std::invalid_argument if min_hit_count == 0.
     */
    void query(
        TermPostings<PositionT> const& term_postings,
        PositionT                      window_length,
        std::size_t                    min_hit_count,
        std::vector<HitInterval>&      hit_intervals
    ) {
        hit_intervals.clear();

        if( min_hit_count == 0 ) {
            throw std::invalid_argument( "HitCollector::query: min_hit_count must be >= 1" );
        }
        if( term_postings.list_count() == 0 || min_hit_count > term_postings.list_count() ) {
            return;
        }
        if( window_length >= static_cast<PositionT>( RingCap ) ) {
            throw std::invalid_argument(
                "HitCollector::query: window_length >= RingCap; increase the RingCap template parameter"
            );
        }
        if( term_postings.list_count() > kMaxLists ) {
            throw std::invalid_argument(
                "HitCollector::query: list_count() > kMaxLists; ListBitset supports at most "
                + std::to_string( kMaxLists ) + " lists"
            );
        }

        // Build a fresh tournament tree, sized to the actual number of input lists, and
        // seed it with one bin per posting list.
        genesis::util::container::TournamentTree<PositionT> tree( term_postings.list_count() );
        tree.build( term_postings.lists() );

        // Ring buffer state
        // ring_head and ring_tail are monotonically increasing;
        // access via `& (RingCap-1)`, stored in ring_mod.
        // Invariant: ring_tail - ring_head == num entries currently in the window.
        std::size_t ring_head = 0;
        std::size_t ring_tail = 0;
        std::size_t const ring_mod = RingCap - 1;

        // Ring access helpers
        auto ring_size  = [&]() noexcept
        {
            return ring_tail - ring_head;
        };
        auto ring_front = [&]() noexcept -> WindowEntry&
        {
            return ring_[ ring_head & ring_mod ];
        };
        auto ring_back = [&]() noexcept -> WindowEntry&
        {
            return ring_[ (ring_tail - 1) & ring_mod ];
        };

        // Hit tracking state
        ListBitset  list_hits{};
        bool        in_valid     = false;
        PositionT   region_start = {};
        PositionT   region_end   = {};
        std::size_t peak_count   = 0;

        // Streaming k-way merge
        while( !tree.empty() ) {

            // Read the globally minimum value, then advance past it
            auto const& top = tree.top();
            PositionT const   val      = top.value;
            std::size_t const list_idx = top.list_index;
            tree.pop();

            // Add to ring buffer (tail-check dedup)
            // Equal values from different lists arrive consecutively from the min-heap,
            // so a matching entry is always at the tail of the ring.
            if( ring_size() > 0 && ring_back().value == val ) {
                ring_back().bits.set( list_idx );
            } else {
                assert( ring_size() < RingCap );
                ring_[ ring_tail & ring_mod ] = { val, {} };
                ring_[ ring_tail & ring_mod ].bits.set( list_idx );
                ++ring_tail;
            }

            // Shrink from front: evict values outside [val-window_length, val]
            // Subtraction is safe: heap yields sorted values, so val >= ring_front().value.
            while( ring_size() > 0 && val - ring_front().value > window_length ) {
                ++ring_head;
            }

            // Recompute window join from scratch
            // At most window_length+1 entries; each join is 4 uint64_t operations.
            list_hits.reset();
            for( std::size_t j = ring_head; j < ring_tail; ++j ) {
                list_hits.join_with( ring_[ j & ring_mod ].bits );
            }

            // Update maximal valid-region tracking
            int const hits = list_hits.hit_count();
            if( hits >= static_cast<int>( min_hit_count ) ) {
                if( !in_valid ) {
                    in_valid     = true;
                    region_start = ring_front().value;
                    peak_count   = 0;
                }
                region_end = ring_back().value;
                if( static_cast<std::size_t>( hits ) > peak_count ) {
                    peak_count = static_cast<std::size_t>( hits );
                }
            } else if( in_valid ) {
                in_valid = false;
                hit_intervals.push_back({ region_start, region_end, peak_count });
            }
        }

        // Flush any open valid region at end of stream
        if( in_valid ) {
            hit_intervals.push_back({ region_start, region_end, peak_count });
        }

        // Sort by peak_count descending: strongest regions first
        std::sort( hit_intervals.begin(), hit_intervals.end(), []( HitInterval const& a, HitInterval const& b ) {
            return a.peak_count > b.peak_count;
        });
    }

    /**
     * @brief Find all maximal hit intervals; returns results by value.
     *
     * Prefer the out-parameter overload to reuse the result vector across calls.
     */
    std::vector<HitInterval> query(
        TermPostings<PositionT> const& term_postings,
        PositionT                      window_length,
        std::size_t                    min_hit_count
    ) {
        std::vector<HitInterval> hit_intervals;
        query( term_postings, window_length, min_hit_count, hit_intervals );
        return hit_intervals;
    }

    // -------------------------------------------------------------------------
    //     Private types and members
    // -------------------------------------------------------------------------

private:

    /**
     * @brief Fixed-capacity bitset tracking which posting lists are represented in the window.
     *
     * One bit per posting list; capacity is kMaxLists (must be a multiple of 64). Stored as
     * kWords = kMaxLists/64 uint64_t words. All operations are branchless word-level instructions
     * unrolled at compile time via index_sequence fold expressions — updating kMaxLists
     * automatically scales all operations with no manual changes required here.
     */
    struct ListBitset
    {
        static constexpr std::size_t kWords = kMaxLists / 64;
        std::array<std::uint64_t, kWords> w = {};

        void set( std::size_t i ) noexcept
        {
            w[i >> 6] |= std::uint64_t{1} << (i & 63u);
        }

        void join_with( ListBitset const& o ) noexcept
        {
            // Fold over indices 0..kWords-1 at compile time; expands to kWords OR assignments
            // with no runtime loop. Equivalent to w[0] |= o.w[0]; w[1] |= o.w[1]; ... for any kWords.
            [&]<std::size_t... I>( std::index_sequence<I...> ) noexcept {
                (( w[I] |= o.w[I] ), ...);
            }( std::make_index_sequence<kWords>{} );
        }

        void reset() noexcept
        {
            w = {};
        }

        int hit_count() const noexcept
        {
            // Fold over indices 0..kWords-1 at compile time; expands to a sum of kWords popcount
            // calls with no runtime loop. The + fold collapses to a single integer at compile time.
            return [&]<std::size_t... I>( std::index_sequence<I...> ) noexcept {
                return static_cast<int>(( std::popcount( w[I] ) + ... ));
            }( std::make_index_sequence<kWords>{} );
        }
    };

    struct WindowEntry
    {
        PositionT  value;
        ListBitset bits;
    };

    std::array<WindowEntry, RingCap> ring_{};
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_HIT_COLLECTOR_H_
