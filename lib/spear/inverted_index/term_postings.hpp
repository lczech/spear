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
#include <string>
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
        std::size_t found       = 0;
        std::size_t empty       = 0;
        std::size_t hard_capped = 0;
        std::size_t soft_capped = 0;
    };

    // -------------------------------------------------------------------------
    //     Construction
    // -------------------------------------------------------------------------

    /**
     * @brief Construct with a pool of @p max_terms decode buffers.
     *
     * @p max_terms must be at least the total number of add() calls that will be made
     * before the next clear(). kEmpty, kHardCapped, and kSoftCapped results do not consume a
     * permanent slot but they do require a temporary slot during decoding, so size this to
     * the total number of add() calls (typically S, the number of queried terms).
     *
     * @param max_posting_list_length  Forwarded to InvertedIndex::postings() on every add();
     *        terms whose posting count exceeds this are treated as kSoftCapped, skipping
     *        decode entirely. Default 0 means no soft cap. See set_max_posting_list_length().
     */
    explicit TermPostings( std::size_t max_terms, std::uint64_t max_posting_list_length = 0 )
        : lists_( max_terms )
        , max_posting_list_length_( max_posting_list_length )
    {
        deferred_.reserve( max_terms );
    }

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    /**
     * @brief Set the soft cap on posting list length, forwarded to InvertedIndex::postings().
     *
     * Terms whose posting count exceeds this are treated as kSoftCapped, skipping decode
     * entirely. A value of 0 disables the soft cap. Unlike the index's own (hard) cap, this
     * can be changed at any time without rebuilding the index, e.g. to experiment with
     * different limits across runs.
     */
    void set_max_posting_list_length( std::uint64_t max_posting_list_length ) noexcept
    {
        max_posting_list_length_ = max_posting_list_length;
    }

    /**
     * @brief Current soft cap on posting list length; 0 means unlimited.
     */
    std::uint64_t max_posting_list_length() const noexcept
    {
        return max_posting_list_length_;
    }

    // -------------------------------------------------------------------------
    //     Accumulation interface
    // -------------------------------------------------------------------------

    /**
     * @brief Reset for a new query round, keeping buffer allocations intact.
     *
     * Throws std::runtime_error if add_deferred() was called but flush_deferred() was not,
     * since that indicates deferred terms would be silently lost.
     */
    void clear()
    {
        if( !deferred_.empty() ) {
            throw std::runtime_error(
                "TermPostings::clear: " + std::to_string( deferred_.size() ) +
                " deferred term(s) were never flushed; call flush_deferred() first"
            );
        }
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
     * kEmpty, kHardCapped, and kSoftCapped results are counted in the histogram but do not
     * claim a slot; the same slot is reused by the next add() call.
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
        // Decode into the next slot. For kEmpty/kHardCapped/kSoftCapped the slot is cleared
        // by postings() but count_ is not incremented, so the slot is reused by the next
        // add() call.
        auto const status = index.postings( term_idx, lists_[count_], max_posting_list_length_ );

        switch( status ) {
            case PostingsStatus::kFound:
                ++count_;
                ++hist_.found;
                break;
            case PostingsStatus::kEmpty:
                ++hist_.empty;
                break;
            case PostingsStatus::kHardCapped:
                ++hist_.hard_capped;
                break;
            case PostingsStatus::kSoftCapped:
                ++hist_.soft_capped;
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
    //     Deferred (Prefetch-Pipelined) Accumulation
    // -------------------------------------------------------------------------

    /**
     * @brief Enqueue @p term_idx for later decoding and prefetch its offset-table entry.
     *
     * Stores the term index in an internal buffer and fires a non-blocking prefetch hint
     * for the offset-table entry so it arrives in L2 cache before flush_deferred() runs.
     * Combine with flush_deferred() for a prefetch-pipelined two-phase decode:
     *
     * @code
     *   tp.clear();
     *   for( auto term : query_terms ) {
     *       tp.add_deferred( index, term );   // called e.g. during k-mer extraction
     *   }
     *   tp.flush_deferred( index );           // blob prefetch + decode in one pipelined pass
     * @endcode
     *
     * May be freely mixed with direct add() calls in the same query round.
     */
    void add_deferred(
        InvertedIndex<PositionT> const& index,
        term_index_type                 term_idx
    ) {
        deferred_.push_back( term_idx );
        index.prefetch_offset( term_idx );
    }

    /**
     * @brief Decode all terms enqueued via add_deferred(), using a prefetch-pipelined loop.
     *
     * Runs two interleaved passes over the deferred term list:
     *   - A priming loop fires blob prefetches for the first @p prefetch_distance terms.
     *   - A combined loop decodes term[i] while prefetching the blob for term[i + prefetch_distance].
     *
     * At steady state, exactly @p prefetch_distance blob-start prefetches are in flight,
     * staying within hardware line-fill buffer capacity.  Each blob prefetch covers 256 bytes
     * (four cache lines), which handles most short posting lists outright and primes the
     * hardware sequential prefetcher for longer ones.
     *
     * Clears the deferred queue on completion.  @p prefetch_distance = 8 is a good default
     * for typical DRAM latency (~150 ns) and short posting list decode times (~10–30 ns each);
     * tune upward if lists are longer and decode is slower.
     */
    void flush_deferred(
        InvertedIndex<PositionT> const& index,
        std::size_t                     prefetch_distance = 8
    ) {
        auto const n     = deferred_.size();
        auto const prime = std::min( prefetch_distance, n );

        // Prime: fire blob prefetches for the first prefetch_distance terms
        for( std::size_t i = 0; i < prime; ++i ) {
            index.prefetch_blob( deferred_[i] );
        }

        // Pipelined loop: decode term[i], prefetch blob for term[i + prefetch_distance]
        for( std::size_t i = 0; i < n; ++i ) {
            if( i + prefetch_distance < n ) {
                index.prefetch_blob( deferred_[i + prefetch_distance] );
            }
            add( index, deferred_[i] );
        }

        deferred_.clear();
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

    /**
     * @brief The successfully decoded (kFound) posting lists, in add() order.
     *
     * Only the first list_count() entries of the underlying pool are valid; slots beyond
     * that may hold stale data from a previous, larger round. The returned span covers
     * exactly the valid entries.
     */
    std::span<std::vector<PositionT> const> lists() const noexcept
    {
        return std::span<std::vector<PositionT> const>( lists_.data(), count_ );
    }

    // -------------------------------------------------------------------------
    //     Data members
    // -------------------------------------------------------------------------

private:

    std::vector<std::vector<PositionT>> lists_;
    std::vector<term_index_type>        deferred_;
    std::size_t                         count_ = 0;
    Stats                               hist_  = {};
    std::uint64_t                       max_posting_list_length_ = 0;
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_TERM_POSTINGS_H_
