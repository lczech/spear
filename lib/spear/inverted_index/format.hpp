#ifndef SPEAR_LIB_INVERTED_INDEX_FORMAT_H_
#define SPEAR_LIB_INVERTED_INDEX_FORMAT_H_

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
 * @brief On-disk format constants and layout structs for the inverted index file format.
 *
 * @file
 * @ingroup util
 */

#include <cassert>
#include <cstdint>

namespace spear::inverted_index {

// =================================================================================================
//     Inverted Index File Format
// =================================================================================================

/**
 * @brief Magic number identifying an inverted index file.
 *
 * Encodes the ASCII string "SPEARII\0" as a little-endian uint64.
 * A reader should reject any file whose footer magic does not match this value.
 */
static constexpr std::uint64_t INVERTED_INDEX_MAGIC = 0x0049495241455053ULL;

/**
 * @brief File format version.
 *
 * Incremented whenever the on-disk layout changes in a backwards-incompatible way.
 * Readers must reject files whose footer version does not equal this constant,
 * unless they implement the necessary backwards compatibility.
 */
static constexpr std::uint64_t INVERTED_INDEX_VERSION = 1;

/**
 * @brief Fixed-size footer written at the very end of every inverted index file.
 *
 * A reader bootstraps by seeking to `file_size - sizeof(InvertedIndexFooter)`, reading
 * this struct, validating magic and version, and then using the offsets to locate the
 * offset table and posting-list blobs.
 *
 * File layout (all multi-byte integers are native-endian):
 * @code
 *   [compressed posting blobs  — variable length, concatenated]
 *   [offset table header       — 3 × uint64: size, width_a, width_b]
 *   [offset table storage      — BitpackedPairVector raw words]
 *   [footer                    — sizeof(InvertedIndexFooter) bytes]
 * @endcode
 *
 * The offset table stores (byte_offset, count) pairs for every term index.
 * An empty term has count == 0; a capped term has count == max_postings_per_term + 1.
 *
 * All fields are uint64_t for alignment and simplicity.
 */
struct InvertedIndexFooter
{
    /// File format magic number. Must equal INVERTED_INDEX_MAGIC.
    std::uint64_t magic;

    /// File format version. Must equal INVERTED_INDEX_VERSION; reject the file otherwise.
    std::uint64_t version;

    /// Number of term indices in this index (i.e., size of the offset table).
    std::uint64_t num_terms;

    /// Byte offset from the start of the file to the first byte of the offset table.
    std::uint64_t offset_table_offset;

    /// Size in bytes of the entire offset table (2-word header + raw bit-packed storage).
    std::uint64_t offset_table_bytes;

    /// Cap threshold used when building. Terms with more than this many unique postings
    /// were discarded; they appear in the offset table with count == max_postings_per_term + 1.
    std::uint64_t max_postings_per_term;

    /// Bit width of stored positions: either 32 or 64. Determines which PFor codec
    /// (pfor_decode_delta1<uint32_t> vs pfor_decode_delta1<uint64_t>) the reader must use.
    std::uint64_t position_bits;

    /// Highest position value that was ever passed to InvertedIndexBuilder::add(),
    /// regardless of whether it was ultimately stored (e.g., for capped terms).
    std::uint64_t max_position;
};

static_assert( sizeof(InvertedIndexFooter) == 64, "InvertedIndexFooter must be exactly 64 bytes" );

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_FORMAT_H_
