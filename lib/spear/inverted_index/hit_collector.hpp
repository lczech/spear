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

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Hit Collector
// =================================================================================================

/**
 * @brief Finds maximal hit intervals across a TermPostings collection.
 *
 * Given a TermPostings collection, a window length @p L, and a minimum hit count @p M,
 * finds every maximal interval [v_left, v_right] in posting-value space such that a
 * sliding window of size L can be placed anywhere within it and still see at least M
 * distinct posting lists represented.
 *
 * Internally runs a streaming k-way merge over the S sorted posting lists via a min-heap,
 * maintaining a small ring buffer (size @p RingCap, a power of two, such that RingCap > L)
 * of (value, Bitset256) window entries. All scratch state (heap, ring buffer) is held as
 * members and reused across calls to avoid per-query allocations.
 *
 * @tparam PositionT  Unsigned integer type for posting positions.
 * @tparam RingCap    Ring buffer capacity; must be a power of two and satisfy RingCap > L
 *                    for every L passed to query(). Default 16 supports L up to 15.
 *                    Use 128 for L up to 127, etc.
 */
template<typename PositionT = std::uint64_t, std::size_t RingCap = 16>
class HitCollector
{
    static_assert( (RingCap & (RingCap - 1)) == 0, "RingCap must be a power of two" );

public:

    // -------------------------------------------------------------------------
    //     Typedefs
    // -------------------------------------------------------------------------

    using position_type = PositionT;

    /**
     * @brief A maximal interval [v_left, v_right] where the hit condition is continuously met.
     *
     * Every window of size L placed within this interval captures at least M distinct lists.
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
     * @brief Find all maximal hit intervals; results are appended to @p out (cleared first),
     * and sorted by peak_count descending, then v_left ascending.
     *
     * @param coll  Populated TermPostings collection.
     * @param L     Window length: a window [x, x+L] must cover values from >= M lists.
     * @param M     Minimum number of distinct lists required; must be >= 1.
     * @param out   Output vector; cleared then filled with results. Sorted by peak hit count.
     *
     * Returns immediately with an empty @p out if M > coll.list_count().
     * Throws std::invalid_argument if M == 0.
     */
    void query(
        TermPostings<PositionT> const& coll,
        PositionT                      L,
        std::size_t                    M,
        std::vector<HitInterval>&      out
    ) {
        out.clear();

        if( M == 0 ) {
            throw std::invalid_argument( "HitCollector::query: M must be >= 1" );
        }
        if( coll.list_count() == 0 || M > coll.list_count() ) {
            return;
        }
        if( L >= static_cast<PositionT>( RingCap ) ) {
            throw std::invalid_argument(
                "HitCollector::query: L >= RingCap; increase the RingCap template parameter"
            );
        }
        if( coll.list_count() > 256 ) {
            throw std::invalid_argument(
                "HitCollector::query: list_count() > 256; Bitset256 supports at most 256 lists"
            );
        }

        // Seed the min-heap (one entry per non-empty posting list)
        heap_.clear();
        heap_.reserve( coll.list_count() );
        for( std::size_t i = 0; i < coll.list_count(); ++i ) {
            auto const sp = coll.list_at( i );
            if( !sp.empty() ) {
                heap_.push_back({
                    sp[0],
                    sp.data(),
                    sp.data() + sp.size(),
                    i
                });
            }
        }
        std::make_heap( heap_.begin(), heap_.end(), HeapCmp{} );

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
        Bitset256   list_hits{};
        bool        in_valid     = false;
        PositionT   region_start = {};
        PositionT   region_end   = {};
        std::size_t peak_count   = 0;

        // Streaming k-way merge
        while( !heap_.empty() ) {

            // Pop the globally minimum value
            std::pop_heap( heap_.begin(), heap_.end(), HeapCmp{} );
            HeapEntry& e = heap_.back();

            PositionT const   val      = e.value;
            std::size_t const list_idx = e.list_idx;

            // Advance this list's cursor; re-insert if not exhausted, else remove
            ++e.cur;
            if( e.cur != e.end ) {
                e.value = *e.cur;
                std::push_heap( heap_.begin(), heap_.end(), HeapCmp{} );
            } else {
                heap_.pop_back();
            }

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

            // Shrink from front: evict values outside [val-L, val]
            // Subtraction is safe: heap yields sorted values, so val >= ring_front().value.
            while( ring_size() > 0 && val - ring_front().value > L ) {
                ++ring_head;
            }

            // Recompute window join from scratch
            // At most L+1 entries; each join is 4 uint64_t operations.
            list_hits.reset();
            for( std::size_t j = ring_head; j < ring_tail; ++j ) {
                list_hits.join_with( ring_[ j & ring_mod ].bits );
            }

            // Update maximal valid-region tracking
            int const hits = list_hits.hit_count();
            if( hits >= static_cast<int>( M ) ) {
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
                out.push_back({ region_start, region_end, peak_count });
            }
        }

        // Flush any open valid region at end of stream
        if( in_valid ) {
            out.push_back({ region_start, region_end, peak_count });
        }

        // Sort by peak_count descending: strongest regions first
        std::sort( out.begin(), out.end(), []( HitInterval const& a, HitInterval const& b ) {
            return a.peak_count > b.peak_count;
        });
    }

    /**
     * @brief Find all maximal hit intervals; returns results by value.
     *
     * Prefer the out-parameter overload to reuse the result vector across calls.
     */
    std::vector<HitInterval> query(
        TermPostings<PositionT> const& coll,
        PositionT                      L,
        std::size_t                    M
    ) {
        std::vector<HitInterval> out;
        query( coll, L, M, out );
        return out;
    }

    // -------------------------------------------------------------------------
    //     Private types and members
    // -------------------------------------------------------------------------

private:

    /**
     * @brief Fixed 256-bit bitset for tracking up to 256 distinct term list memberships.
     *
     * Stored as four `uint64_t` words. All operations are branchless word-level instructions.
     *
     * To extend beyond 256 lists, four things need to change:
     *  1. Grow `w`: increase the `std::array` size (N words → 64*N bit capacity).
     *  2. `join_with`: add one `w[k] |= o.w[k]` line per new word.
     *  3. `popcount`: add one `std::popcount(w[k])` term per new word.
     *  4. The guard in `HitCollector::query`: update the `list_count() > 256` threshold
     *     to match the new capacity (64 * array size).
     * `set()` already uses `i >> 6` / `i & 63u` generically and needs no change.
     * `HeapEntry::list_idx` is `std::size_t` and also needs no change.
     */
    struct Bitset256
    {
        std::array<std::uint64_t, 4> w = {};

        void set( std::size_t i ) noexcept
        {
            w[i >> 6] |= std::uint64_t{1} << (i & 63u);
        }

        void join_with( Bitset256 const& o ) noexcept
        {
            w[0] |= o.w[0];
            w[1] |= o.w[1];
            w[2] |= o.w[2];
            w[3] |= o.w[3];
        }

        void reset() noexcept
        {
            w = {};
        }

        int hit_count() const noexcept
        {
            return
                std::popcount( w[0] ) + std::popcount( w[1] ) +
                std::popcount( w[2] ) + std::popcount( w[3] )
            ;
        }
    };

    struct WindowEntry
    {
        PositionT value;
        Bitset256 bits;
    };

    struct HeapEntry
    {
        PositionT          value;    // current value at cursor (cached to avoid indirection)
        PositionT const*   cur;      // pointer to current element in the posting list
        PositionT const*   end;      // one past the last element
        std::size_t        list_idx; // which slot in TermPostings
    };

    struct HeapCmp
    {
        bool operator()( HeapEntry const& a, HeapEntry const& b ) const noexcept
        {
            return a.value > b.value; // min-heap: smaller value has higher priority
        }
    };

    std::array<WindowEntry, RingCap> ring_{};
    std::vector<HeapEntry>           heap_;
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_HIT_COLLECTOR_H_
