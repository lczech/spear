#ifndef SPEAR_LIB_INVERTED_INDEX_PFOR_H_
#define SPEAR_LIB_INVERTED_INDEX_PFOR_H_

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

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace spear::inverted_index {

// =================================================================================================
//     32-bit
// =================================================================================================

/**
 * @brief Returns the maximum encoded output buffer size in bytes for @p n uint32_t integers.
 */
[[nodiscard]] std::size_t pfor_bound_32( std::size_t n ) noexcept;

// -----------------------------------------------------------------------------
// Plain P4, uint32_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n arbitrary uint32_t integers from @p in into @p out.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n uint32_t integers from @p in into @p out.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
);

// -----------------------------------------------------------------------------
//     Delta P4, uint32_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n non-decreasing uint32_t integers using delta coding.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_delta_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n non-decreasing uint32_t integers from delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_delta_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
);

// -----------------------------------------------------------------------------
//     Delta1 P4, uint32_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n strictly increasing uint32_t integers using delta-1 coding.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_delta1_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n strictly increasing uint32_t integers from delta-1-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_delta1_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
);

// -----------------------------------------------------------------------------
//     Zigzag P4, uint32_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n uint32_t integers using zigzag-of-delta coding.
 *
 * Best suited for sequences where consecutive elements differ by small signed amounts
 * (i.e., values can go up or down by small steps). Each element is reinterpreted as
 * int32_t, the signed delta to the previous element is computed, and that delta is
 * zigzag-mapped to a small unsigned value before PFor encoding.
 *
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_zigzag_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n uint32_t integers from zigzag-of-delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_zigzag_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
);

// =================================================================================================
//     64-bit
// =================================================================================================

/**
 * @brief Returns the maximum encoded output buffer size in bytes for @p n uint64_t integers.
 */
[[nodiscard]] std::size_t pfor_bound_64( std::size_t n ) noexcept;

// -----------------------------------------------------------------------------
//     Plain P4, uint64_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n arbitrary uint64_t integers from @p in into @p out.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n uint64_t integers from @p in into @p out.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
);

// -----------------------------------------------------------------------------
//     Delta P4, uint64_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n non-decreasing uint64_t integers using delta coding.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_delta_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n non-decreasing uint64_t integers from delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_delta_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
);

// -----------------------------------------------------------------------------
//     Delta1 P4, uint64_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n strictly increasing uint64_t integers using delta-1 coding.
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_delta1_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n strictly increasing uint64_t integers from delta-1-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_delta1_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
);

// -----------------------------------------------------------------------------
//     Zigzag P4, uint64_t
// -----------------------------------------------------------------------------

/**
 * @brief PFor-encode @p n uint64_t integers using zigzag-of-delta coding.
 *
 * Best suited for sequences where consecutive elements differ by small signed amounts
 * (i.e., values can go up or down by small steps). Each element is reinterpreted as
 * int64_t, the signed delta to the previous element is computed, and that delta is
 * zigzag-mapped to a small unsigned value before PFor encoding.
 *
 * @return Number of bytes written to @p out.
 */
[[nodiscard]] std::size_t pfor_encode_zigzag_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
);

/**
 * @brief PFor-decode @p n uint64_t integers from zigzag-of-delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
[[nodiscard]] std::size_t pfor_decode_zigzag_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
);

// =================================================================================================
//     Type-dispatching wrappers (uint32_t / uint64_t)
// =================================================================================================

/**
 * @brief Returns the maximum encoded output buffer size in bytes for @p n integers of type @p T.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_bound( std::size_t n ) noexcept
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_bound requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_bound_32(n);
    } else {
        return pfor_bound_64(n);
    }
}

/**
 * @brief PFor-encode @p n arbitrary integers of type @p T from @p in into @p out.
 * @return Number of bytes written to @p out.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_encode( T* in, std::size_t n, std::uint8_t* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_encode requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_encode_32(in, n, out);
    } else {
        return pfor_encode_64(in, n, out);
    }
}

/**
 * @brief PFor-decode @p n integers of type @p T from @p in into @p out.
 * @return Number of bytes consumed from @p in.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_decode( std::uint8_t* in, std::size_t n, T* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_decode requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_decode_32(in, n, out);
    } else {
        return pfor_decode_64(in, n, out);
    }
}

/**
 * @brief PFor-encode @p n non-decreasing integers of type @p T using delta coding.
 * @return Number of bytes written to @p out.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_encode_delta( T* in, std::size_t n, std::uint8_t* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_encode_delta requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_encode_delta_32(in, n, out);
    } else {
        return pfor_encode_delta_64(in, n, out);
    }
}

/**
 * @brief PFor-decode @p n non-decreasing integers of type @p T from delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_decode_delta( std::uint8_t* in, std::size_t n, T* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_decode_delta requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_decode_delta_32(in, n, out);
    } else {
        return pfor_decode_delta_64(in, n, out);
    }
}

/**
 * @brief PFor-encode @p n strictly increasing integers of type @p T using delta-1 coding.
 * @return Number of bytes written to @p out.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_encode_delta1( T* in, std::size_t n, std::uint8_t* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_encode_delta1 requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_encode_delta1_32(in, n, out);
    } else {
        return pfor_encode_delta1_64(in, n, out);
    }
}

/**
 * @brief PFor-decode @p n strictly increasing integers of type @p T from delta-1-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_decode_delta1( std::uint8_t* in, std::size_t n, T* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_decode_delta1 requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_decode_delta1_32(in, n, out);
    } else {
        return pfor_decode_delta1_64(in, n, out);
    }
}

/**
 * @brief PFor-encode @p n integers of type @p T using zigzag-of-delta coding.
 *
 * Best suited for sequences where consecutive elements differ by small signed amounts.
 * See @ref pfor_encode_zigzag_32 / @ref pfor_encode_zigzag_64 for details.
 *
 * @return Number of bytes written to @p out.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_encode_zigzag( T* in, std::size_t n, std::uint8_t* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_encode_zigzag requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_encode_zigzag_32(in, n, out);
    } else {
        return pfor_encode_zigzag_64(in, n, out);
    }
}

/**
 * @brief PFor-decode @p n integers of type @p T from zigzag-of-delta-coded @p in.
 * @return Number of bytes consumed from @p in.
 */
template<typename T>
[[nodiscard]] inline std::size_t pfor_decode_zigzag( std::uint8_t* in, std::size_t n, T* out )
{
    static_assert(
        std::is_same_v<T, std::uint32_t> || std::is_same_v<T, std::uint64_t>,
        "pfor_decode_zigzag requires uint32_t or uint64_t"
    );
    if constexpr (std::is_same_v<T, std::uint32_t>) {
        return pfor_decode_zigzag_32(in, n, out);
    } else {
        return pfor_decode_zigzag_64(in, n, out);
    }
}

} // namespace spear::inverted_index

#endif // GENESIS_UTIL_BIT_PFOR_H_
