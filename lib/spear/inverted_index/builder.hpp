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
#include <bit>
#include <cassert>
#include <cerrno>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "spear/inverted_index/format.hpp"
#include "spear/inverted_index/pfor.hpp"
#include "spear/inverted_index/stats.hpp"
#include "genesis/util/container/bitpacked_pair_vector.hpp"
#include "genesis/util/container/bitpacked_pair_vector_io.hpp"
#include "genesis/util/io/file_handle.hpp"
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
 *     has been capped, via a sentinel value.
 *
 * Semantics:
 *   - The per-entry counter stores the number of integers currently present in
 *     the compressed PFOR data.
 *   - Pending positions are not counted until flush().
 *   - Once capped, no postings are stored anymore for that term, only its counts.
 *
 * Thread safety:
 *   - add() is thread-safe.
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
        , pending_positions_(get_pending_storage_size_(num_term_indices, cfg.pending_capacity))
        , entries_(num_term_indices)
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
            cfg_.max_postings_per_term != 0 &&
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

    InvertedIndexBuilder( InvertedIndexBuilder&& other ) noexcept
        : cfg_(std::move(other.cfg_))
        , pending_positions_(std::move(other.pending_positions_))
        , entries_(std::move(other.entries_))
        , written_count_(other.written_count_)
        , guards_(std::move(other.guards_))
    {
        // Leave other in a defined state: empty entries, no pointers to free.
        other.written_count_ = 0;
    }

    InvertedIndexBuilder& operator= ( InvertedIndexBuilder const& ) = delete;

    InvertedIndexBuilder& operator= ( InvertedIndexBuilder&& other ) noexcept
    {
        if( this != &other ) {
            // Free raw pointers for our own active entries before overwriting.
            for( std::size_t i = written_count_; i < entries_.size(); ++i ) {
                delete[] entries_[i].active.compressed_positions;
            }
            cfg_               = std::move(other.cfg_);
            pending_positions_ = std::move(other.pending_positions_);
            entries_           = std::move(other.entries_);
            guards_            = std::move(other.guards_);
            written_count_     = other.written_count_;
            other.written_count_ = 0;
        }
        return *this;
    }

    ~InvertedIndexBuilder()
    {
        // Entries [0, written_count_) are in written state: their raw pointers have already
        // been freed by write(). Only free the remaining active entries.
        // Under normal operation, written_count_ should be zero until write() and then jump to
        // entries_.size() at the end of write(), so this loop should usually run either zero
        // times or all times. Other states only happen in case of an exception during write().
        for( std::size_t i = written_count_; i < entries_.size(); ++i ) {
            delete[] entries_[i].active.compressed_positions;
        }
    }

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
            entry.active.pending_count > 0 &&
            pending[static_cast<std::size_t>(entry.active.pending_count) - 1] == position
        ) {
            return;
        }

        // If there is still room in the pending buffer, just append.
        if( entry.active.pending_count < static_cast<pending_count_type>(cfg_.pending_capacity) ) {
            pending[entry.active.pending_count] = position;
            ++entry.active.pending_count;
            return;
        }

        // Otherwise, the new element would overflow the pending buffer, so flush the
        // existing pending values *and* this new value together.
        flush_pending_locked_(index, &position);
    }

    /**
     * @brief Flush all pending buffers into compressed postings.
     *
     * Must not race with add(). That is, finalize() must only be called
     * after all producer threads have finished, and only from one thread.
     */
    void finalize()
    {
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            // auto lock = guards_.get_lock_guard(i);
            if( !is_capped_(entries_[i]) && entries_[i].active.pending_count > 0 ) {
                flush_pending_locked_(i);
            }
        }
    }

    /**
     * @brief Write the index to disk and return statistics.
     *
     * Processes all entries (flushing any remaining pending data), writes compressed posting
     * lists as one contiguous byte blob, then appends a bit-packed offset table and a fixed-size
     * footer. After this call the builder is in the same state as if freshly constructed with the
     * same num_term_indices and Config: all compressed data has been freed.
     *
     * File layout (all multi-byte integers are native-endian):
     * @code
     *   [compressed posting blobs  — variable, concatenated]
     *   [offset table header       — 2 × uint64: width_a, width_b]
     *   [offset table storage      — BitpackedPairVector raw words]
     *   [footer                    — 6 × uint64 = 48 bytes]
     * @endcode
     *
     * The footer is an InvertedIndexFooter (see format.hpp): magic, version, num_terms,
     * offset_table_offset, offset_table_bytes, max_postings_per_term, position_bits.
     * The offset table stores (byte_offset, count) pairs indexed by term index.
     * For capped terms, count == max_postings_per_term + 1; for empty terms, count == 0.
     *
     * @param path Output file path. Throws if the file already exists (unless
     *             genesis::util::core::Options::get().allow_file_overwriting() is set).
     * @return Statistics about the written index.
     */
    InvertedIndexStats write( std::string const& path )
    {
        // We delegate, for easy try-catch wrapping.
        assert( written_count_ == 0 );
        try {
            return write_impl_( path );
        } catch(...) {
            reset_after_failed_write_();
            throw;
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
        return static_cast<std::size_t>(entry.active.compressed_count);
    }

    /**
     * @brief Get the postings for a given term index, decompressed into a vector of positions.
     */
    [[nodiscard]] std::vector<position_type> postings(term_index_type term_index) const
    {
        std::vector<position_type> result;
        postings( term_index, result );
        return result;
    }

    /**
     * @brief Get the postings for a given term index, decompressed into @p target_buf.
     *
     * Resizes @p target_buf to the posting count and fills it with decoded positions.
     * For capped or empty terms @p target_buf is cleared. Reusing @p target_buf across calls
     * avoids repeated allocations on the hot path.
     */
    void postings(term_index_type term_index, std::vector<position_type>& target_buf) const
    {
        auto const index = checked_index_(term_index);
        auto lock = guards_.get_lock_guard(index);

        // Check if the term is capped or has no data.
        Entry const& entry = entries_[index];
        if( is_capped_(entry) || entry.active.compressed_count == 0 ) {
            target_buf.clear();
            return;
        }

        // Decode the compressed postings into buf.
        // Auto-switch between 32 and 64 bit PFOR based on the position_type.
        target_buf.resize( static_cast<std::size_t>( entry.active.compressed_count ));
        pfor_decode_delta1<position_type>(
            entry.active.compressed_positions,
            static_cast<std::size_t>(entry.active.compressed_count),
            target_buf.data()
        );
    }

private:

    // -------------------------------------------------------------------------
    //     Internal structs
    // -------------------------------------------------------------------------

    // Each entry slot is a union of two states:
    //   active:  the normal build state — a raw pointer to PFOR-compressed data plus counters.
    //   written: set by write() once an entry is processed — the byte offset and post count
    //            for the offset table, stored in the same 16 bytes without extra allocation.
    //
    // written_count_ is the boundary: entries [0, written_count_) are in written state,
    // entries [written_count_, entries_.size()) are in active state. The destructor and
    // move operators use this to free only the active entries' raw pointers.
    union Entry
    {
        struct Active
        {
            // Raw pointer to PFOR-compressed posting data. Owned by InvertedIndexBuilder;
            // freed explicitly in cap_entry_, reset_entry_, write(), and the destructor.
            std::uint8_t* compressed_positions;

            // Number of compressed postings in compressed_positions.
            stored_count_type compressed_count;

            // Number of pending positions currently in the pending buffer for this term index.
            pending_count_type pending_count;

            // compressed_positions grows through kBinSizes[]; compressed_bin is the index
            // into that table giving the allocated capacity. Meaningful only when
            // compressed_positions != nullptr. Using uint8_t keeps Entry at 16 bytes (with
            // padding due to pointer alignment) for a wider range of template combinations.
            std::uint8_t compressed_bin;
        } active;

        struct Written
        {
            std::uint64_t offset; // byte offset of this entry's blob in the output file
            std::uint64_t count;  // posting count, or capped sentinel, or 0 for empty
        } written;

        Entry() : active{ nullptr, 0, 0, 0 } {}
        ~Entry() {} // InvertedIndexBuilder owns the raw pointer lifecycle
    };

    // -------------------------------------------------------------------------
    //     Write to file
    // -------------------------------------------------------------------------

    // Body of write(). Separated so write() can wrap it in a clean try/catch.
    InvertedIndexStats write_impl_( std::string const& path )
    {
        // Get a binary output FILE* with safety checks and RAII management.
        auto fp_guard = genesis::util::io::file_output_file(path);
        FILE* const fp = fp_guard.get();

        // Helper lambda to check fwrite() return values and throw on error with a nice message.
        auto checked_fwrite = [&]( void const* ptr, std::size_t bytes ) {
            if( std::fwrite( ptr, 1, bytes, fp ) != bytes ) {
                throw std::runtime_error(
                    "Write error on '" + path + "': " + std::string( std::strerror( errno ))
                );
            }
        };

        // Prepare stats. Histogram indices: 0 = empty, 1..max = count, max+1 = capped.
        // Grown on demand so no-cap mode reflects actual distribution without pre-sizing.
        InvertedIndexStats stats;
        auto hist_inc = [&]( std::uint64_t val ) {
            if( val >= stats.posting_count_histogram.size() ) {
                stats.posting_count_histogram.resize( val + 1, 0 );
            }
            ++stats.posting_count_histogram[val];
        };

        // Reusable encode buffer (heap-allocated, grown as needed).
        std::vector<std::uint8_t> encode_buf;

        // Track the maximum post count seen (needed for width_b when there is no cap).
        std::uint64_t max_post_count = 0;

        // Main blob write loop.
        // As each entry is processed, its raw pointer is freed and the slot is repurposed
        // in-place to hold the (offset, count) pair. written_count_ is updated immediately
        // so the destructor can free remaining active pointers on any exception path.
        std::uint64_t current_offset = 0;
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            Entry& entry = entries_[i];
            std::uint64_t const entry_offset = current_offset;
            std::uint64_t post_count = 0;

            // Capped entries do not have a blob; their pointer is already null from cap_entry_.
            if( is_capped_(entry) ) {
                post_count = static_cast<std::uint64_t>(capped_sentinel());
            } else {
                // Decode existing + merge pending into thread_local buffer.
                auto const& buf = flush_to_buffer_(i, nullptr);
                std::size_t const count = buf.size();

                if( count == 0 ) {
                    // Empty entry: nothing to write, post_count stays 0.
                } else if(
                    cfg_.max_postings_per_term != 0 &&
                    count > static_cast<std::size_t>(cfg_.max_postings_per_term)
                ) {
                    // Pending data pushed an uncapped entry over the cap threshold.
                    post_count = static_cast<std::uint64_t>(capped_sentinel());
                } else {
                    // Regular entry with data: encode and write.
                    std::size_t const bound = pfor_bound<position_type>(count);
                    encode_buf.resize(bound);
                    std::size_t const bytes_written = pfor_encode_delta1<position_type>(
                        buf.data(), count, encode_buf.data()
                    );
                    checked_fwrite( encode_buf.data(), bytes_written );
                    stats.total_blob_bytes += bytes_written;
                    current_offset += bytes_written;
                    post_count = static_cast<std::uint64_t>(count);
                }

                // Free the raw pointer for this entry's active state.
                delete[] entry.active.compressed_positions;
            }

            hist_inc( post_count );
            max_post_count = std::max( max_post_count, post_count );

            // Repurpose the slot in-place: raw pointer is gone, store (offset, count) here.
            // Update written_count_ before the next iteration so the destructor sees
            // a consistent boundary between written and active entries at all times.
            entry.written = { entry_offset, post_count };
            written_count_ = i + 1;
        }

        // Offset table.
        // Determine bit widths from observed maxima. width_a covers byte offsets
        // (max = total blob bytes); width_b covers counts including the capped sentinel.
        auto needed_bits = []( std::uint64_t val ) -> std::size_t {
            return val == 0 ? 1u : static_cast<std::size_t>(std::bit_width(val));
        };
        std::size_t const width_a = needed_bits( current_offset );
        std::size_t const width_b = ( cfg_.max_postings_per_term != 0 )
            ? needed_bits( static_cast<std::uint64_t>(capped_sentinel()))
            : needed_bits( max_post_count )
        ;
        if( width_a + width_b > 64 ) {
            throw std::overflow_error(
                "Offset table bit widths (" + std::to_string(width_a) + " + " +
                std::to_string(width_b) + ") exceed 64 bits"
            );
        }

        // Construct the offset table directly from the repurposed entry slots.
        using BPV = genesis::util::container::BitpackedPairVector<
            std::uint64_t, std::uint64_t, std::uint64_t
        >;
        BPV offset_table( entries_.size(), width_a, width_b );
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            offset_table.set_unchecked( i, entries_[i].written.offset, entries_[i].written.count );
        }

        // Compute offset table byte size: 3-word header + raw storage words.
        std::size_t const num_storage_words = offset_table.storage().data().size();
        std::uint64_t const offset_table_offset = current_offset;
        std::uint64_t const offset_table_bytes  =
            3 * sizeof(std::uint64_t) +                   // header: size, width_a, width_b
            num_storage_words * sizeof(BPV::storage_type) // raw bitpacked storage
        ;

        // Write offset table (header + raw storage) via genesis I/O helper; verify byte count.
        using namespace genesis::util::container;
        [[maybe_unused]] std::size_t const written = write_bitpacked_pair_vector(
            offset_table, fp
        );
        assert( written == static_cast<std::size_t>( offset_table_bytes ));

        // Write footer with metadata about the index.
        InvertedIndexFooter const footer{
            INVERTED_INDEX_MAGIC,
            INVERTED_INDEX_VERSION,
            static_cast<std::uint64_t>(entries_.size()),
            offset_table_offset,
            offset_table_bytes,
            static_cast<std::uint64_t>(cfg_.max_postings_per_term),
            static_cast<std::uint64_t>(sizeof(position_type) * 8)
        };
        checked_fwrite( &footer, sizeof(footer) );

        // Reset all entries to active-empty state, restoring the builder for reuse.
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            entries_[i].active = { nullptr, 0, 0, 0 };
        }
        written_count_ = 0;

        return stats;
    }

    // Called by write() on any exception: frees remaining active entries and resets all
    // slots to active-empty state, leaving the builder in the same state as after construction.
    void reset_after_failed_write_() noexcept
    {
        // Entries [written_count_, N) are still active — free their raw pointers.
        for( std::size_t i = written_count_; i < entries_.size(); ++i ) {
            delete[] entries_[i].active.compressed_positions;
        }
        // Reset all slots (including already-written ones) to active-empty state.
        for( std::size_t i = 0; i < entries_.size(); ++i ) {
            entries_[i].active = { nullptr, 0, 0, 0 };
        }
        written_count_ = 0;
    }

    // -------------------------------------------------------------------------
    //     Entry processing
    // -------------------------------------------------------------------------

    [[nodiscard]] inline std::size_t checked_index_(term_index_type term_index) const
    {
        auto const index = static_cast<std::size_t>(term_index);
        if (index >= entries_.size()) {
            throw std::out_of_range("term index out of range");
        }
        return index;
    }

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
        return
            cfg_.max_postings_per_term != 0 &&
            entry.active.compressed_count == capped_sentinel()
        ;
    }

    void cap_entry_(Entry& entry) noexcept
    {
        delete[] entry.active.compressed_positions;
        entry.active.compressed_positions = nullptr;
        entry.active.compressed_bin       = 0;
        entry.active.compressed_count     = capped_sentinel();
        entry.active.pending_count        = 0;
    }

    /**
     * @brief Flush one entry's pending data into its compressed storage.
     *
     * Precondition: caller holds the per-term guard.
     */
    void flush_pending_locked_(std::size_t index)
    {
        flush_pending_locked_(index, nullptr);
    }

    /**
     * @brief Flush one entry with an optional extra position into its compressed storage.
     *
     * Precondition: exclusive access to this entry is guaranteed — either by holding the
     * per-term guard (add() path) or because no concurrent writers exist (finalize() path).
     */
    void flush_pending_locked_(std::size_t index, position_type const* extra_position)
    {
        Entry& entry = entries_[index];
        assert(!is_capped_(entry));

        // Nothing to do.
        if( entry.active.pending_count == 0 && extra_position == nullptr ) {
            return;
        }

        // Get the merged existing + pending + extra data into a thread_local buffer.
        auto const& decode_buf = flush_to_buffer_(index, extra_position);
        std::size_t const count = decode_buf.size();

        // Cap if we now exceed the max postings per term.
        if(
            cfg_.max_postings_per_term != 0 &&
            count > static_cast<std::size_t>(cfg_.max_postings_per_term)
        ) {
            cap_entry_(entry);
            return;
        }
        if( count > static_cast<std::size_t>(std::numeric_limits<stored_count_type>::max()) ) {
            throw std::overflow_error("posting count exceeds StoredCountT");
        }

        // Recompress and store in the entry.
        // Encode into a thread_local staging buffer sized for the worst case, then memcpy
        // only the bytes actually written into the entry's buffer. The entry buffer grows
        // geometrically by sqrt(2) when needed, so most flushes are a plain memcpy with no
        // malloc/free — O(log(final_compressed_size)) allocations per term over its lifetime.
        static thread_local std::vector<std::uint8_t> encode_buf;
        encode_buf.resize( pfor_bound<position_type>( count ));
        std::size_t const bytes_written = pfor_encode_delta1<position_type>(
            decode_buf.data(), count, encode_buf.data()
        );

        // Check if new allocation is needed or current bin is sufficient.
        if(
            !entry.active.compressed_positions ||
            bytes_written > kBinSizes[entry.active.compressed_bin]
        ) {
            std::uint8_t const new_bin = size_to_bin_( bytes_written );
            delete[] entry.active.compressed_positions;
            entry.active.compressed_positions = new std::uint8_t[kBinSizes[new_bin]];
            entry.active.compressed_bin       = new_bin;
        }

        // Store the compressed data and count in the entry.
        std::memcpy( entry.active.compressed_positions, encode_buf.data(), bytes_written );
        entry.active.compressed_count = static_cast<stored_count_type>(count);
    }

    /**
     * @brief Merge existing compressed data and pending buffer into a thread_local buffer.
     *
     * Decodes compressed_positions (if any), appends pending positions and an optional extra
     * value, then sorts the new tail, merges with the sorted existing prefix, and deduplicates.
     * Resets entry.pending_count to 0. Returns a reference to the (thread_local) result buffer.
     *
     * The reference is valid until the next call to flush_to_buffer_ on any entry in this thread.
     * Precondition: entry is not capped.
     */
    [[nodiscard]] std::vector<position_type> const& flush_to_buffer_(
        std::size_t index,
        position_type const* extra_position
    ) {
        Entry& entry = entries_[index];
        assert( !is_capped_(entry) );
        assert( entry.active.pending_count <= static_cast<pending_count_type>(cfg_.pending_capacity) );
        assert( entry.active.compressed_count == 0 || entry.active.compressed_positions != nullptr );

        // Get all involved counts of entries.
        std::size_t const existing_count = static_cast<std::size_t>(entry.active.compressed_count);
        std::size_t const pending_count  = static_cast<std::size_t>(entry.active.pending_count);
        std::size_t const extra_count    = extra_position ? 1u : 0u;
        std::size_t const total_count    = existing_count + pending_count + extra_count;

        // Use a local buffer to avoid reallocations, amortized over time.
        static thread_local std::vector<position_type> buffer;
        // buffer.clear();
        buffer.resize(total_count);

        // Decode old compressed postings into the front.
        if( existing_count > 0 ) {
            (void) pfor_decode_delta1<position_type>(
                entry.active.compressed_positions,
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
        entry.active.pending_count = 0;

        // Sort only the newly appended tail; the existing head is already sorted.
        auto const middle = buffer.begin() + static_cast<std::ptrdiff_t>(existing_count);
        std::sort(middle, buffer.end());

        // Merge old sorted prefix with new sorted suffix, and remove duplicates.
        std::inplace_merge(buffer.begin(), middle, buffer.end());
        buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());

        return buffer;
    }

    // -------------------------------------------------------------------------
    //     Bin Table for Compressed Buffer Sizes
    // -------------------------------------------------------------------------

    // Allocation sizes for compressed_positions follow a sqrt(2) progression, giving at most
    // ~1.41x over-allocation instead of the ~2x waste of pure power-of-two rounding.
    // Even-indexed bins are exact powers of two; odd-indexed bins are ceil(2^(i/2)),
    // computed via the integer approximation 181/128 ≈ sqrt(2).
    // 80 bins cover up to ~777 GB, sufficient for uint32_t StoredCountT + uint64_t positions
    // with a conservative 2x headroom. The static_assert below enforces this at compile time.

    static constexpr std::size_t kNumBins = 80;

    static constexpr std::array<std::size_t, kNumBins> kBinSizes = []() {
        std::array<std::size_t, kNumBins> b{};
        for( std::size_t i = 0; i < kNumBins; ++i ) {
            if( i % 2 == 0 ) {
                b[i] = std::size_t{1} << (i / 2);
            } else {
                // ceil(2^(i/2) * sqrt(2)); 181/128 ≈ 1.4140625 ≈ sqrt(2), +127 for ceiling
                b[i] = (std::size_t{181} * (std::size_t{1} << (i / 2)) + 127) / 128;
            }
        }
        return b;
    }();

    // Verify the table is large enough for the worst-case compressed size.
    // Conservative bound: PFOR-delta on sorted input never expands beyond 2x uncompressed.
    static constexpr std::size_t kMaxPostingsForType_ =
        static_cast<std::size_t>(std::numeric_limits<stored_count_type>::max()) - 1;
    static constexpr std::size_t kMaxCompressedBound_ =
        kMaxPostingsForType_ * sizeof(position_type) * 2;
    static_assert(
        kBinSizes.back() >= kMaxCompressedBound_,
        "kBinSizes table too small for this StoredCountT / position_type combination; "
        "increase kNumBins"
    );

    // Return the index of the smallest bin whose size is >= bytes.
    [[nodiscard]] static std::uint8_t size_to_bin_( std::size_t bytes )
    {
        auto const it = std::lower_bound( kBinSizes.begin(), kBinSizes.end(), bytes );
        assert( it != kBinSizes.end() ); // static_assert above should prevent this
        return static_cast<std::uint8_t>( it - kBinSizes.begin() );
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    Config cfg_;

    // Contiguous storage for all pending buffers,
    // indexed by [term_index * pending_capacity + offset].
    std::vector<position_type> pending_positions_;

    // Term index of compressed postings. Each slot is a union of an active entry
    // (raw pointer + counters) and a written entry (offset + count); see Entry above.
    std::vector<Entry> entries_;

    // Number of leading entries already converted to written state by write().
    // Entries [0, written_count_) hold (offset, count) pairs; their raw pointers are freed.
    // Entries [written_count_, entries_.size()) are in active state and may have raw pointers.
    // The destructor and move operators use this boundary for correct pointer cleanup.
    std::size_t written_count_ = 0;

    // Thread-safe access for concurrent insertion into entries.
    mutable genesis::util::threading::ConcurrentVectorGuard guards_;
};

} // namespace spear::inverted_index

#endif // include guard
