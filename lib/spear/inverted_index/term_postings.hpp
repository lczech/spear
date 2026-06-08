#ifndef SPEAR_LIB_INVERTED_INDEX_TERM_POSTINGS_H_
#define SPEAR_LIB_INVERTED_INDEX_TERM_POSTINGS_H_

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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Term Postings
// =================================================================================================

/**
 * @brief Accumulates decoded posting lists for a set of query terms.
 *
 * Owns a fixed pool of decode buffers (one per expected term). Buffers are reused across
 * queries via clear(), which resets the count but keeps the vector capacities allocated.
 * Terms with empty or capped posting lists are skipped (not assigned a buffer slot) but
 * counted in the histogram for diagnostics.
 *
 * Typical usage:
 * @code
 *     // num_terms is the max number of terms that will be queried together
 *     TermPostings<uint64_t> tp(num_terms);
 *
 *     // For each query:
 *     tp.clear();
 *     for (auto term_idx : query_terms) {
 *         tp.add(inverted_index, term_idx);
 *     }
 *     // pass tp to HitCollector::query()
 *     // examine tp.stats() for diagnostics
 * @endcode
 *
 * @tparam PositionT  Unsigned integer type for posting positions; must match the InvertedIndex.
 */
template<typename PositionT = std::uint64_t>
class TermPostings
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs
    // -------------------------------------------------------------------------

    using position_type   = PositionT;
    using term_index_type = typename InvertedIndex<PositionT>::term_index_type;
    using PostingsStatus  = typename InvertedIndex<PositionT>::PostingsStatus;

    // -------------------------------------------------------------------------
    //     Stats
    // -------------------------------------------------------------------------

    /**
     * @brief Counts of each PostingsStatus outcome across all add() calls.
     */
    struct Stats
    {
        std::size_t found  = 0;
        std::size_t empty  = 0;
        std::size_t capped = 0;
    };

    // -------------------------------------------------------------------------
    //     Construction
    // -------------------------------------------------------------------------

    /**
     * @brief Construct with a pool of @p max_terms decode buffers.
     *
     * @p max_terms must be at least the total number of add() calls that will be made
     * before the next clear(). kEmpty and kCapped results do not consume a permanent slot
     * but they do require a temporary slot during decoding, so size this to the total
     * number of add() calls (typically S, the number of queried terms).
     */
    explicit TermPostings( std::size_t max_terms )
        : lists_( max_terms )
    {}

    // -------------------------------------------------------------------------
    //     Accumulation interface
    // -------------------------------------------------------------------------

    /**
     * @brief Reset for a new query round, keeping buffer allocations intact.
     */
    void clear()
    {
        for( std::size_t i = 0; i < count_; ++i ) {
            lists_[i].clear();
        }
        count_ = 0;
        hist_  = {};
    }

    /**
     * @brief Decode the posting list for @p term_idx and add it to the collection.
     *
     * Throws std::runtime_error if the pool is full (i.e. add() has been called more
     * times than max_terms allows). kFound results claim a slot and increment list_count().
     * kEmpty and kCapped results are counted in the histogram but do not claim a slot;
     * the same slot is reused by the next add() call.
     *
     * Returns the PostingsStatus for caller inspection.
     */
    PostingsStatus add(
        InvertedIndex<PositionT> const& index,
        term_index_type                 term_idx
    ) {
        if( count_ >= lists_.size() ) {
            throw std::runtime_error(
                "TermPostings::add: pool capacity exceeded; increase max_terms"
            );
        }
        // Decode into the next slot. For kEmpty/kCapped the slot is cleared by postings()
        // but count_ is not incremented, so the slot is reused by the next add() call.
        auto const status = index.postings( term_idx, lists_[count_] );

        switch( status ) {
            case PostingsStatus::kFound:
                ++count_;
                ++hist_.found;
                break;
            case PostingsStatus::kEmpty:
                ++hist_.empty;
                break;
            case PostingsStatus::kCapped:
                ++hist_.capped;
                break;
        }
        return status;
    }

    /**
     * @brief Add a pre-sorted posting list directly, without going through an InvertedIndex.
     *
     * The list is moved into the next pool slot; the caller's vector is left empty.
     * An empty vector is recorded as kEmpty and does not consume a pool slot.
     * Throws std::runtime_error if the pool is full and @p list is non-empty.
     *
     * The list must be sorted in ascending order; no validation is performed.
     */
    PostingsStatus add( std::vector<PositionT> list )
    {
        if( list.empty() ) {
            ++hist_.empty;
            return PostingsStatus::kEmpty;
        }
        if( count_ >= lists_.size() ) {
            throw std::runtime_error(
                "TermPostings::add: pool capacity exceeded; increase max_terms"
            );
        }
        lists_[count_] = std::move( list );
        ++count_;
        ++hist_.found;
        return PostingsStatus::kFound;
    }

    /**
     * @brief Add a pre-sorted posting list from a span, copying the values.
     *
     * An empty span is recorded as kEmpty and does not consume a pool slot.
     * Throws std::runtime_error if the pool is full and @p list is non-empty.
     *
     * The list must be sorted in ascending order; no validation is performed.
     */
    PostingsStatus add( std::span<PositionT const> list )
    {
        if( list.empty() ) {
            ++hist_.empty;
            return PostingsStatus::kEmpty;
        }
        if( count_ >= lists_.size() ) {
            throw std::runtime_error(
                "TermPostings::add: pool capacity exceeded; increase max_terms"
            );
        }
        lists_[count_].assign( list.begin(), list.end() );
        ++count_;
        ++hist_.found;
        return PostingsStatus::kFound;
    }

    // -------------------------------------------------------------------------
    //     Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Number of successfully decoded (kFound) posting lists.
     */
    std::size_t list_count() const noexcept
    {
        return count_;
    }

    /**
     * @brief Posting list for the i-th successfully decoded term (zero-indexed).
     */
    std::span<const PositionT> list_at( std::size_t i ) const noexcept
    {
        assert( i < count_ );
        return lists_[i];
    }

    /**
     * @brief Stats for PostingsStatus outcomes across all add() calls since last clear().
     */
    Stats const& stats() const noexcept
    {
        return hist_;
    }

    // -------------------------------------------------------------------------
    //     Data members
    // -------------------------------------------------------------------------

private:

    std::vector<std::vector<PositionT>> lists_;
    std::size_t                         count_ = 0;
    Stats                               hist_  = {};
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_TERM_POSTINGS_H_
