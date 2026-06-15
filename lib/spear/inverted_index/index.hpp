#ifndef SPEAR_LIB_INVERTED_INDEX_INDEX_H_
#define SPEAR_LIB_INVERTED_INDEX_INDEX_H_

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
 * @brief Reader and accessor for the on-disk inverted index file format.
 *
 * @file
 * @ingroup util
 */

#include "spear/inverted_index/format.hpp"
#include "spear/inverted_index/pfor.hpp"
#include "genesis/util/container/bitpacked_pair_vector.hpp"
#include "genesis/util/container/bitpacked_pair_vector_io.hpp"
#include "genesis/util/io/file_handle.hpp"

#include <cassert>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Inverted Index Reader / Accessor
// =================================================================================================

/**
 * @brief Inverted index, to look up posting lists for indexed terms.
 *
 * Opens an index file produced by InvertedIndexBuilder and provides
 * random-access lookup of posting lists by term index.
 *
 * The offset table is always loaded into RAM at open time. The posting blob (compressed posting
 * lists per term) is handled according to the chosen OpenMode:
 *
 *   - **kLoadAll**: the entire blob is read into a `vector<uint8_t>` at open time.
 *     All postings() calls are pointer arithmetic + PFor decode; no further I/O.
 *     Recommended for cluster nodes with sufficient RAM.
 *
 *   - **kPread**: only the footer and offset table are loaded; the blob stays on disk.
 *     Each postings() call issues a `pread()` syscall for exactly the compressed bytes
 *     of the requested term. Thread-safe. Suitable when the blob does not fit in RAM,
 *     and the I/O overhead is acceptable (e.g., on a fast SSD or NVMe).
 *
 * Both modes are safe to call from multiple threads concurrently: in kLoadAll mode all
 * state is read-only after open(); in kPread mode each thread uses a thread_local decode
 * buffer and `pread()` is inherently thread-safe.
 *
 * @tparam PositionT
 *     Unsigned integer type for decoded posting positions. Must be `uint32_t` or `uint64_t`,
 *     matching the type used when the file was built. The footer's `position_bits` field is
 *     checked against `sizeof(PositionT) * 8` at open time; a mismatch throws.
 */
template<typename PositionT = std::uint64_t>
class InvertedIndex
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Enums
    // -------------------------------------------------------------------------

    using position_type  = PositionT;
    using term_index_type = std::size_t;

    static_assert(
        std::is_same_v<PositionT, std::uint32_t> || std::is_same_v<PositionT, std::uint64_t>,
        "InvertedIndex requires PositionT = uint32_t or uint64_t"
    );

    // -------------------------------------------------------------------------
    //     Open Modes and Status
    // -------------------------------------------------------------------------

    /**
     * @brief Controls how the posting blob is accessed after open().
     */
    enum class OpenMode
    {
        /// Load the entire blob into RAM at open time. Zero I/O per postings() call.
        kLoadAll,

        /// Keep blob on disk; issue one pread() per postings() call.
        kPread
    };

    /**
     * @brief Result status returned by postings().
     */
    enum class PostingsStatus : std::uint8_t
    {
        /// The term has no postings; result is empty.
        kEmpty,

        /// The posting list was decoded into the output buffer.
        kFound,

        /// The term exceeded max_postings_per_term and was discarded at build time;
        /// result is empty. This is a "hard cap" baked into the index itself.
        kHardCapped,

        /// The term's posting count exceeded the caller-supplied max_posting_list_length
        /// for this call; result is empty, but the data is still present in the index
        /// (unlike kHardCapped). This is a "soft cap" applied per query.
        kSoftCapped
    };

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    InvertedIndex() = default;

    ~InvertedIndex()
    {
        // Close the file descriptor kept open for pread mode
        if( fd_ >= 0 ) {
            ::close( fd_ );
        }
    }

    // Copying is disabled: the blob can be hundreds of GB and fd cannot be trivially duplicated.
    InvertedIndex( InvertedIndex const& ) = delete;
    InvertedIndex& operator=( InvertedIndex const& ) = delete;

    InvertedIndex( InvertedIndex&& other ) noexcept
    {
        swap_( other );
    }

    InvertedIndex& operator=( InvertedIndex&& other ) noexcept
    {
        InvertedIndex tmp( std::move( other ));
        swap_( tmp );
        return *this;
    }

    // -------------------------------------------------------------------------
    //     Open
    // -------------------------------------------------------------------------

    /**
     * @brief Open an inverted index file.
     *
     * This loads the footer, offset table, and (optionally) the blob.  May be called again to
     * reopen a different file; any previously held blob or file descriptor is released first.
     *
     * @param path     Path to the index file produced by InvertedIndexBuilder::write().
     * @param mode     kLoadAll (default) or kPread.
     * @throws std::runtime_error on any I/O error or format mismatch.
     */
    void open( std::string const& path, OpenMode mode = OpenMode::kLoadAll )
    {
        // Reset any previously opened state
        if( fd_ >= 0 ) {
            ::close( fd_ );
            fd_ = -1;
        }
        blob_.clear();
        blob_.shrink_to_fit();
        offset_table_ = OffsetTable{};
        footer_ = InvertedIndexFooter{};
        mode_ = mode;

        // Open file (throws on failure, RAII-managed)
        auto fp_guard = genesis::util::io::file_input_file( path );
        FILE* const fp = fp_guard.get();

        // Seek to footer position and read the fixed-size footer struct
        if( fseeko( fp, -static_cast<off_t>( sizeof( InvertedIndexFooter )), SEEK_END ) != 0 ) {
            throw std::runtime_error( "Cannot seek to footer in '" + path + "'" );
        }
        if( std::fread( &footer_, sizeof( footer_ ), 1, fp ) != 1 ) {
            throw std::runtime_error( "Cannot read footer from '" + path + "'" );
        }

        // Validate magic number and format version
        if( footer_.magic != INVERTED_INDEX_MAGIC ) {
            throw std::runtime_error(
                "'" + path + "' is not an inverted index file (bad magic)"
            );
        }
        if( footer_.version != INVERTED_INDEX_VERSION ) {
            throw std::runtime_error(
                "'" + path + "' has format version " + std::to_string( footer_.version ) +
                " but this reader expects version " + std::to_string( INVERTED_INDEX_VERSION )
            );
        }

        // Validate position type: file must match the template parameter
        constexpr std::uint64_t expected_bits = static_cast<std::uint64_t>( sizeof( PositionT ) * 8 );
        if( footer_.position_bits != expected_bits ) {
            throw std::runtime_error(
                "'" + path + "' was built with " + std::to_string( footer_.position_bits ) +
                "-bit positions but InvertedIndex<" +
                ( sizeof( PositionT ) == 4 ? "uint32_t" : "uint64_t" ) +
                "> expects " + std::to_string( expected_bits ) + "-bit positions"
            );
        }

        // Seek to offset table and load it into the in-memory BPV
        if( fseeko( fp, static_cast<off_t>( footer_.offset_table_offset ), SEEK_SET ) != 0 ) {
            throw std::runtime_error( "Cannot seek to offset table in '" + path + "'" );
        }
        offset_table_ = genesis::util::container::read_bitpacked_pair_vector<
            std::uint64_t, std::uint64_t, std::uint64_t
        >( fp );
        if( offset_table_.size() != static_cast<std::size_t>( footer_.num_terms )) {
            throw std::runtime_error(
                "'" + path + "': offset table size " + std::to_string( offset_table_.size() ) +
                " does not match footer num_terms " + std::to_string( footer_.num_terms )
            );
        }

        if( mode_ == OpenMode::kLoadAll ) {
            // Load the entire posting blob into memory for zero-I/O query access
            blob_.resize( static_cast<std::size_t>( footer_.offset_table_offset ));
            if( !blob_.empty() ) {
                if( fseeko( fp, 0, SEEK_SET ) != 0 ) {
                    throw std::runtime_error( "Cannot seek to blob start in '" + path + "'" );
                }
                if( std::fread( blob_.data(), 1, blob_.size(), fp ) != blob_.size() ) {
                    throw std::runtime_error( "Cannot read posting blob from '" + path + "'" );
                }
            }
            // fp_guard closes the file here
        } else {
            // pread mode: open a raw fd for thread-safe positional reads; fp_guard closes FILE*
            fd_ = ::open( path.c_str(), O_RDONLY );
            if( fd_ < 0 ) {
                throw std::runtime_error(
                    "Cannot open '" + path + "' for pread: " + std::string( std::strerror( errno ))
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    //     Posting Queries
    // -------------------------------------------------------------------------

    /**
     * @brief Decode and return the posting list for @p term_index.
     *
     * Returns the PostingsStatus as well to indicate the status of the result.
     * Thread-safe in both OpenMode variants.
     */
    [[nodiscard]] std::pair<std::vector<PositionT>, PostingsStatus>
    postings( term_index_type term_index ) const
    {
        std::vector<PositionT> buf;
        PostingsStatus const status = postings( term_index, buf );
        return { std::move( buf ), status };
    }

    /**
     * @brief Decode the posting list for @p term_index into @p target_buf.
     *
     * Resizes @p target_buf to the posting count and fills it with decoded positions.
     * For empty, capped, or soft-capped terms @p target_buf is cleared.  Reusing
     * @p target_buf across calls avoids repeated allocations on the hot path.
     *
     * @param max_posting_list_length  If non-zero, terms whose posting count exceeds this
     *        value are treated as kSoftCapped: the decode is skipped entirely (the count
     *        is available from the offset table alone), but unlike kHardCapped, the data
     *        remains in the index and would be returned with a higher or zero limit.
     *        This allows experimenting with tighter caps at query time without rebuilding
     *        the index. Default 0 means no soft cap is applied.
     *
     * Returns PostingsStatus::kFound, kEmpty, kHardCapped, or kSoftCapped accordingly.
     * Thread-safe in both OpenMode variants.
     */
    PostingsStatus
    postings(
        term_index_type term_index, std::vector<PositionT>& target_buf,
        std::uint64_t max_posting_list_length = 0
    ) const {
        assert( footer_.magic == INVERTED_INDEX_MAGIC );
        assert( footer_.num_terms == offset_table_.size() );

        // Validate term index
        if( term_index >= static_cast<term_index_type>( footer_.num_terms )) {
            throw std::out_of_range(
                "term index " + std::to_string( term_index ) + " out of range "
                "(num_terms=" + std::to_string( footer_.num_terms ) + ")"
            );
        }

        // Look up the byte offset and count for this term
        auto const [byte_offset, count] = offset_table_[term_index];

        // Empty and hard-capped terms carry no blob data
        if( count == 0 ) {
            target_buf.clear();
            return PostingsStatus::kEmpty;
        }
        if( footer_.max_postings_per_term != 0 && count >= capped_sentinel() ) {
            target_buf.clear();
            return PostingsStatus::kHardCapped;
        }
        if( max_posting_list_length != 0 && count > max_posting_list_length ) {
            target_buf.clear();
            return PostingsStatus::kSoftCapped;
        }

        // Resize the output buffer to hold exactly the decoded positions
        std::size_t const n = static_cast<std::size_t>( count );
        target_buf.resize( n );

        if( mode_ == OpenMode::kLoadAll ) {
            // Decode directly from the in-memory blob via pointer arithmetic
            assert( byte_offset < blob_.size() );
            [[maybe_unused]] std::size_t const decoded_bytes = pfor_decode_delta1<PositionT>(
                blob_.data() + byte_offset, n, target_buf.data()
            );
            assert( [&]{
                std::uint64_t const next_blob_offset =
                    ( term_index + 1 < footer_.num_terms )
                    ? offset_table_.first_at( term_index + 1 )
                    : footer_.offset_table_offset
                ;
                return decoded_bytes == static_cast<std::size_t>( next_blob_offset - byte_offset );
            }() );
        } else {
            // Derive compressed byte length from consecutive offsets
            std::uint64_t const next_offset =
                ( term_index + 1 < footer_.num_terms )
                ? offset_table_.first_at( term_index + 1 )
                : footer_.offset_table_offset;
            std::size_t const blob_bytes = static_cast<std::size_t>( next_offset - byte_offset );
            assert( blob_bytes > 0 );

            // Read compressed bytes into a thread-local buffer to avoid per-call allocation
            static thread_local std::vector<std::uint8_t> read_buf;
            read_buf.resize( blob_bytes );
            ssize_t const n_read = ::pread(
                fd_, read_buf.data(), blob_bytes, static_cast<off_t>( byte_offset )
            );
            if( n_read != static_cast<ssize_t>( blob_bytes )) {
                throw std::runtime_error(
                    "pread failed (expected " + std::to_string( blob_bytes ) +
                    " bytes): " + std::string( std::strerror( errno ))
                );
            }

            // Decode the compressed posting list into the output buffer
            [[maybe_unused]] std::size_t const decoded_bytes = pfor_decode_delta1<PositionT>(
                read_buf.data(), n, target_buf.data()
            );
            assert( decoded_bytes == blob_bytes );
        }
        return PostingsStatus::kFound;
    }

    // -------------------------------------------------------------------------
    //     Per-Term Metadata
    // -------------------------------------------------------------------------

    /**
     * @brief Return the posting count for @p term_index.
     *
     * Returns 0 for empty terms and capped_sentinel() for capped terms,
     * matching the sentinel convention used by InvertedIndexBuilder.
     */
    [[nodiscard]] std::uint64_t posting_count( term_index_type term_index ) const
    {
        if( term_index >= static_cast<term_index_type>( footer_.num_terms )) {
            throw std::out_of_range( "term index out of range" );
        }
        return offset_table_.second_at( term_index );
    }

    /**
     * @brief Return true if @p term_index was capped during construction.
     */
    [[nodiscard]] bool is_capped( term_index_type term_index ) const
    {
        if( term_index >= static_cast<term_index_type>( footer_.num_terms )) {
            throw std::out_of_range( "term index out of range" );
        }
        return
            footer_.max_postings_per_term != 0 &&
            offset_table_.second_at( term_index ) == capped_sentinel()
        ;
    }

    // -------------------------------------------------------------------------
    //     Prefetch Hints
    // -------------------------------------------------------------------------

    /**
     * @brief Prefetch the offset-table entry for @p term_index into L2 cache.
     *
     * Call this in a tight loop over all query term indices before calling postings()
     * or add_deferred() / flush_deferred() on a TermPostings collection.  By the time
     * the subsequent decode loop starts, the offset-table words should be in cache.
     *
     * No-op on non-GCC/Clang compilers.  Asserts that @p term_index is in range.
     */
    void prefetch_offset( term_index_type term_index ) const noexcept
    {
        assert( term_index < static_cast<term_index_type>( footer_.num_terms ));
        offset_table_.prefetch( term_index );
    }

    /**
     * @brief Prefetch the first 256 bytes of the posting-list blob for @p term_index.
     *
     * Reads the blob byte-offset from the (hopefully already cached) offset-table entry
     * and issues four sequential cache-line prefetch hints.  Covers short posting lists
     * outright; for longer ones the hardware sequential prefetcher takes over after the
     * first few cache lines establish a stride.
     *
     * No-op in kPread mode (blob not in RAM) and on non-GCC/Clang compilers.
     * Asserts that @p term_index is in range.
     */
    void prefetch_blob( term_index_type term_index ) const noexcept
    {
        assert( term_index < static_cast<term_index_type>( footer_.num_terms ));
        if( blob_.empty() ) {
            return;
        }
        auto const base = blob_.data() + offset_table_.first_unchecked( term_index );
        #if defined(__GNUC__) || defined(__clang__)
            __builtin_prefetch( base,       0, 0 );
            __builtin_prefetch( base + 64,  0, 0 );
            __builtin_prefetch( base + 128, 0, 0 );
            __builtin_prefetch( base + 192, 0, 0 );
        #endif
    }

    // -------------------------------------------------------------------------
    //     Index Metadata
    // -------------------------------------------------------------------------

    /**
     * @brief Number of term indices in this index (size of the offset table).
     */
    [[nodiscard]] std::size_t num_terms() const noexcept
    {
        return static_cast<std::size_t>( footer_.num_terms );
    }

    /**
     * @brief Cap threshold used when building the index.
     *
     * Terms with more unique postings than this were discarded.
     */
    [[nodiscard]] std::uint64_t max_postings_per_term() const noexcept
    {
        return footer_.max_postings_per_term;
    }

    /**
     * @brief The sentinel count value stored for capped terms.
     *
     * Equal to max_postings_per_term() + 1. posting_count() returns this value
     * for any term that was capped during construction.
     */
    [[nodiscard]] std::uint64_t capped_sentinel() const noexcept
    {
        return footer_.max_postings_per_term + 1;
    }

    /**
     * @brief Highest position value that was ever passed to InvertedIndexBuilder::add(),
     * regardless of whether it was ultimately stored (e.g., for capped terms).
     */
    [[nodiscard]] std::uint64_t max_position() const noexcept
    {
        return footer_.max_position;
    }

    /**
     * @brief Return true if the index has been opened, i.e., if open() has been called.
     */
    [[nodiscard]] bool is_open() const noexcept
    {
        // Reusing the magic field as a check sentinel here.
        return footer_.magic == INVERTED_INDEX_MAGIC;
    }

    // -------------------------------------------------------------------------
    //     Private Helpers
    // -------------------------------------------------------------------------

private:

    void swap_( InvertedIndex& other ) noexcept
    {
        using std::swap;
        swap( footer_,       other.footer_ );
        swap( offset_table_, other.offset_table_ );
        swap( blob_,         other.blob_ );
        swap( fd_,           other.fd_ );
        swap( mode_,         other.mode_ );
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

    using OffsetTable = genesis::util::container::BitpackedPairVector<
        std::uint64_t, std::uint64_t, std::uint64_t
    >;

    // Footer loaded from disk; magic field doubles as an "is open" sentinel
    InvertedIndexFooter footer_{};

    // Offset table: maps term_index --> ( byte_offset_in_blob, posting_count )
    OffsetTable offset_table_;

    // kLoadAll mode: owns the entire compressed posting blob
    std::vector<std::uint8_t> blob_;

    // pread mode: raw file descriptor for thread-safe positional reads; -1 when unused
    int fd_ = -1;

    OpenMode mode_ = OpenMode::kLoadAll;
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_INDEX_H_
