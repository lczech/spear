#ifndef SPEAR_LIB_INVERTED_INDEX_BUILDER_H_
#define SPEAR_LIB_INVERTED_INDEX_BUILDER_H_

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
 * @brief
 *
 * @file
 * @ingroup util
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "spear/inverted_index/pfor.hpp"
#include "genesis/util/threading/concurrent_vector_guard.hpp"

namespace spear::inverted_index {

// =================================================================================================
//     Inverted Index Builder
// =================================================================================================

/**
 * @brief Dynamic inverted-index builder with PFOR-compressed posting lists.
 *
 * This models a concurrent data structure for efficiently building an inverted index
 * from a stream of (term_index, position) pairs. We do not store document indices here,
 * and assume that those can be retrieved by working backwards from the positions each term
 * appears in, by assuming positions to be from concatenated documents. Term indices also are
 * dense indices into the storage array here.
 *
 * Design:
 *   - For each term index [0, t], keep a small pending buffer of newly inserted positions.
 *   - All pending buffers are stored in one large preallocated array of
 *     num_term_indices * pending_capacity positions.
 *   - When the pending buffer fills, merge it with the already-compressed postings,
 *     deduplicate, and recompress.
 *   - If the number of unique stored postings exceeds the configured maximum,
 *     discard all postings for that term index permanently and remember only that it
 *     has been capped, via a sentiel value.
 *
 * Semantics:
 *   - The per-entry counter stores the number of integers currently present in
 *     the compressed PFOR data.
 *   - Pending positions are not counted until flush().
 *   - Once capped, no postings are stored anymore for that term, only its counts.
 *
 * Thread safety:
 *   - add() and add_many() are thread-safe.
 *   - finalize() must only be called after all producer threads have finished.
 *
 * @tparam PositionT
 *     Type for positions of terms (e.g., positions or regions in a genome) in the reference
 *     documents. Must be 32 or 64 bit, as supported by the selected PFOR codec.
 *     Switching here will automatically switch the codec used for compression as well.
 *
 * @tparam StoredCountT
 *     Integer type used for the per-entry stored-posting counter. This counts
 *     the number of integers currently present in the compressed posting list.
 *     It also stores the capped sentinel value max_postings_per_term + 1.
 *     Default: uint16_t. That is, no more than 2^16 stored postings per term index.
 *
 * @tparam PendingCountT
 *     Integer type used for the size of the pending buffer for one entry.
 *     The configured pending capacity must fit into this type.
 *     Default: uint8_t. That is, no more than 255 pending positions per term index.
 */
template<
    typename PositionT     = std::uint64_t,
    typename StoredCountT  = std::uint16_t,
    typename PendingCountT = std::uint8_t
>
class InvertedIndexBuilder
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    /**
     * @brief Type for indices into the array of terms (e.g., k-mers in dense 2-bit encoding).
     *
     * This is the index into the entry term array. When using, e.g., canonical k-mer encoding,
     * this is the canonical encoding.
     */
    using term_index_type = std::size_t;

    // Aliases for template parameters.
    using position_type      = PositionT;
    using stored_count_type  = StoredCountT;
    using pending_count_type = PendingCountT;

    static_assert(std::is_integral_v<term_index_type>);
    static_assert(std::is_integral_v<position_type>);
    static_assert(std::is_integral_v<stored_count_type>);
    static_assert(std::is_integral_v<pending_count_type>);
    static_assert(
        sizeof(position_type) == 4 || sizeof(position_type) == 8,
        "position_type must be 32 or 64 bit to be able to work with PFOR"
    );

    // -------------------------------------------------------------------------
    //     Public Helper Structs
    // -------------------------------------------------------------------------

    struct Config
    {
        /**
         * @brief Number of pending positions buffered per term (e.g., per k-mer) before flush.
         *
         * Must fit into PendingCountT.
         */
        std::size_t pending_capacity = 16;

        /**
         * @brief Maximum number of unique postings to keep for one term (e.g., per k-mer).
         *
         * Once a flush would produce more than this many unique postings,
         * all postings for that term are discarded permanently.
         */
        stored_count_type max_postings_per_term = 1024;

        /**
         * @brief Number of internal guard buckets for concurrent insertion.
         *
         * If zero, a heuristic based on number of entries and hardware threads
         * is used by ConcurrentVectorGuard.
         */
        std::size_t num_guard_buckets = 0;
    };

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    InvertedIndexBuilder() = default;

    explicit InvertedIndexBuilder( std::size_t num_term_indices, Config cfg = Config{} )
        : cfg_(cfg)
        , entries_(num_term_indices)
        , pending_positions_(get_pending_storage_size_(num_term_indices, cfg.pending_capacity))
        , guards_(
            cfg.num_guard_buckets == 0
            ? genesis::util::threading::ConcurrentVectorGuard(
                    num_term_indices,
                    std::max<std::size_t>(1, std::thread::hardware_concurrency())
                )
            : genesis::util::threading::ConcurrentVectorGuard(cfg.num_guard_buckets)
        )
    {
        if(
            cfg_.pending_capacity >
            static_cast<std::size_t>(std::numeric_limits<pending_count_type>::max())
        ) {
            throw std::invalid_argument(
                "pending_capacity does not fit into PendingCountT"
            );
        }
        if(
            cfg_.max_postings_per_term >=
            static_cast<stored_count_type>(std::numeric_limits<stored_count_type>::max())
        ) {
            throw std::invalid_argument(
                "max_postings_per_term too large for StoredCountT "
                "(needs one extra value for capped sentinel)"
            );
        }
    }

    InvertedIndexBuilder( InvertedIndexBuilder const& ) = delete;
    InvertedIndexBuilder( InvertedIndexBuilder&& )      = default;

    InvertedIndexBuilder& operator= ( InvertedIndexBuilder const& ) = delete;
    InvertedIndexBuilder& operator= ( InvertedIndexBuilder&& )      = default;

    ~InvertedIndexBuilder() = default;

    // -------------------------------------------------------------------------
    //     Modifiers
    // -------------------------------------------------------------------------

    /**
     * @brief Thread-safe insertion of one (term_index, position) pair.
     */
    void add( term_index_type term_index, position_type position )
    {
        // Lock the entry and check its basics. If it is already capped, we can just return.
        std::size_t const index = checked_index_(term_index);
        auto lock = guards_.get_lock_guard(index);
        Entry& entry = entries_[index];
        if( is_capped_(entry) ) {
            return;
        }

        // Check for duplicates in the pending buffer, just by checking the last,
        // as a fast return point. This assumes that different threads do not tend
        // to write to the same indices very often.
        position_type* const pending = pending_begin_(index);
        if(
            entry.pending_count > 0 &&
            pending[static_cast<std::size_t>(entry.pending_count) - 1] == position
        ) {
            return;
        }

        // If there is still room in the pending buffer, just append.
        if( entry.pending_count < static_cast<pending_count_type>(cfg_.pending_capacity) ) {
            pending[entry.pending_count] = position;
            ++entry.pending_count;
            return;
        }

        // Otherwise, the new element would overflow the pending buffer, so flush the
        // existing pending values *and* this new value together.
        flush_locked_(index, &position);
    }

    /**
     * @brief Flush all pending buffers into compressed postings.
     *
     * Must not race with add() / add_many(). That is, finalize() must only be called
     * after all producer threads have finished, and only from one thread.
     */
    void finalize()
    {
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            // auto lock = guards_.get_lock_guard(i);
            if( !is_capped_(entries_[i]) && entries_[i].pending_count > 0 ) {
                flush_locked_(i);
            }
        }
    }

    // -------------------------------------------------------------------------
    //     Properties and Getters
    // -------------------------------------------------------------------------

    /**
     * @brief Get the configuration of this builder.
     */
    [[nodiscard]] Config const& config() const noexcept
    {
        return cfg_;
    }

    /**
     * @brief Get the number of term indices, i.e., the size of the entry array.
     */
    [[nodiscard]] std::size_t num_term_indices() const noexcept
    {
        return entries_.size();
    }

    /**
     * @brief Check if a term index is capped.
     */
    [[nodiscard]] bool is_capped(term_index_type term_index) const
    {
        auto const index = checked_index_(term_index);
        auto lock = guards_.get_lock_guard(index);
        return is_capped_(entries_[index]);
    }

    /**
     * @brief Get the capped sentinel value.
     */
    [[nodiscard]] inline stored_count_type capped_sentinel() const noexcept
    {
        return static_cast<stored_count_type>(cfg_.max_postings_per_term + 1);
    }

    /**
     * @brief Get the number of postings for a given term index.
     */
    [[nodiscard]] std::size_t posting_count(term_index_type term_index) const
    {
        auto const index = checked_index_(term_index);
        auto lock = guards_.get_lock_guard(index);

        Entry const& entry = entries_[index];
        if( is_capped_(entry) ) {
            return capped_sentinel();
        }
        return static_cast<std::size_t>(entry.compressed_count);
    }

    /**
     * @brief Get the postings for a given term index, decompressed into a vector of positions.
     */
    [[nodiscard]] std::vector<position_type> postings(term_index_type term_index) const
    {
        auto const index = checked_index_(term_index);
        auto lock = guards_.get_lock_guard(index);

        // Check if the term is capped or has no data.
        Entry const& entry = entries_[index];
        if( is_capped_(entry) || entry.compressed_count == 0 ) {
            return {};
        }

        // Decode the compressed postings into a vector and return it.
        // Auto-switch between 32 and 64 bit PFOR based on the position_type.
        std::vector<position_type> result(entry.compressed_count);
        (void) pfor_decode_delta1<position_type>(
            entry.compressed_positions.get(),
            static_cast<std::size_t>(entry.compressed_count),
            result.data()
        );
        return result;
    }

    // -------------------------------------------------------------------------
    //     Internal Functions
    // -------------------------------------------------------------------------

private:

    struct Entry
    {
        // Compressed postings for the term index, stored as a byte array from PFor.
        std::unique_ptr<std::uint8_t[]> compressed_positions;

        // Number of compressed postings in the compressed_positions array.
        stored_count_type compressed_count = 0;

        // Number of pending positions currently in the pending buffer for this term index.
        pending_count_type pending_count = 0;
    };

    [[nodiscard]] static std::size_t get_pending_storage_size_(
        std::size_t num_term_indices,
        std::size_t pending_capacity
    ) {
        if (pending_capacity == 0) {
            return 0;
        }
        if (num_term_indices > std::numeric_limits<std::size_t>::max() / pending_capacity) {
            throw std::overflow_error("pending storage size overflows size_t");
        }
        return num_term_indices * pending_capacity;
    }

    [[nodiscard]] inline std::size_t checked_index_(term_index_type term_index) const
    {
        auto const index = static_cast<std::size_t>(term_index);
        if (index >= entries_.size()) {
            throw std::out_of_range("term index out of range");
        }
        return index;
    }

    [[nodiscard]] inline position_type* pending_begin_(std::size_t term_index) noexcept
    {
        return pending_positions_.data() + term_index * cfg_.pending_capacity;
    }

    [[nodiscard]] inline position_type const* pending_begin_(std::size_t term_index) const noexcept
    {
        return pending_positions_.data() + term_index * cfg_.pending_capacity;
    }

    [[nodiscard]] inline bool is_capped_(Entry const& entry) const noexcept
    {
        return entry.compressed_count == capped_sentinel();
    }

    void cap_entry_(Entry& entry) noexcept
    {
        entry.compressed_positions.reset();
        entry.compressed_count = capped_sentinel();
        entry.pending_count = 0;
    }

    /**
     * @brief Flush one entry.
     *
     * Precondition: caller holds the per-term guard.
     */
    void flush_locked_(std::size_t index)
    {
        flush_locked_(index, nullptr);
    }

    /**
     * @brief Flush one entry with an optional extra position.
     *
     * Precondition: caller holds the per-term guard.
     */
    void flush_locked_(std::size_t index, position_type const* extra_position)
    {
        Entry& entry = entries_[index];
        assert(!is_capped_(entry));

        // Nothing to do
        if( entry.pending_count == 0 && extra_position == nullptr ) {
            return;
        }

        // Get all involved counts of entries
        std::size_t const existing_count = static_cast<std::size_t>(entry.compressed_count);
        std::size_t const pending_count  = static_cast<std::size_t>(entry.pending_count);
        std::size_t const extra_count    = extra_position ? 1u : 0u;
        std::size_t const total_count    = existing_count + pending_count + extra_count;

        // Use a local buffer to avoid reallocations, amortized over time.
        static thread_local std::vector<position_type> buffer;
        // buffer.clear();
        buffer.resize(total_count);

        // Decode old compressed postings into the front.
        if( existing_count > 0 ) {
            (void) pfor_decode_delta1<position_type>(
                entry.compressed_positions.get(),
                existing_count,
                buffer.data()
            );
        }

        // Copy pending block directly behind the decoded postings.
        position_type const* const pending = pending_begin_(index);
        if( pending_count > 0 ) {
            std::copy(
                pending,
                pending + pending_count,
                buffer.begin() + static_cast<std::ptrdiff_t>(existing_count)
            );
        }

        // Append optional extra value at the end.
        if( extra_position ) {
            assert( existing_count + pending_count < buffer.size() );
            buffer[existing_count + pending_count] = *extra_position;
        }

        // Pending data is now consumed.
        entry.pending_count = 0;

        // Sort only the newly appended tail; the existing head is already sorted.
        auto const middle = buffer.begin() + static_cast<std::ptrdiff_t>(existing_count);
        std::sort(middle, buffer.end());

        // Merge old sorted prefix with new sorted suffix, and remove duplicates.
        std::inplace_merge(buffer.begin(), middle, buffer.end());
        buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

        // Cap if needed.
        if( buffer.size() > static_cast<std::size_t>(cfg_.max_postings_per_term) ) {
            cap_entry_(entry);
            return;
        }

        // Recompress and store in the entry.
        std::size_t const bound = pfor_bound<position_type>(buffer.size());
        auto new_bytes = std::make_unique<std::uint8_t[]>(bound);
        (void) pfor_encode_delta1<position_type>(
            buffer.data(),
            buffer.size(),
            new_bytes.get()
        );
        if( buffer.size() > static_cast<std::size_t>(std::numeric_limits<stored_count_type>::max()) ) {
            throw std::overflow_error("posting count exceeds StoredCountT");
        }
        entry.compressed_positions = std::move(new_bytes);
        entry.compressed_count = static_cast<stored_count_type>(buffer.size());
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    Config cfg_;

    // Term index of compressed postings
    std::vector<Entry> entries_;

    // Contiguous storage for all pending buffers,
    // indexed by [term_index * pending_capacity + offset].
    std::vector<position_type> pending_positions_;

    mutable genesis::util::threading::ConcurrentVectorGuard guards_;
};

} // namespace spear::inverted_index

#endif // include guard
