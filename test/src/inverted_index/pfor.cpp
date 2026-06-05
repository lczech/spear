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
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/inverted_index/pfor.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

using namespace spear::inverted_index;

// =================================================================================================
//     Helpers
// =================================================================================================

template <typename UInt>
std::vector<UInt> make_patterned_values_( std::size_t n )
{
    std::vector<UInt> v;
    v.reserve( n );

    for( std::size_t i = 0; i < n; ++i ) {
        auto const x = static_cast<UInt>(
            ( i % 7 == 0 ) ? i :
            ( i % 7 == 1 ) ? ( i * 3 ) :
            ( i % 7 == 2 ) ? ( i * 17 ) :
            ( i % 7 == 3 ) ? ( i * i ) :
            ( i % 7 == 4 ) ? ( ( i << 5 ) ^ ( i * 13 ) ) :
            ( i % 7 == 5 ) ? ( ( UInt{1} << ( ( i % ( sizeof(UInt) * 8 - 2 )) + 1 ) ) - 1 ) :
                             ( i * 257 + 11 )
        );
        v.push_back( x );
    }

    return v;
}

template <typename UInt>
std::vector<UInt> make_monotone_values_( std::size_t n )
{
    std::vector<UInt> v;
    v.reserve( n );

    UInt cur = 0;
    for( std::size_t i = 0; i < n; ++i ) {
        UInt const step = static_cast<UInt>(
            ( i % 8 == 0 ) ? 0 :
            ( i % 8 == 1 ) ? 1 :
            ( i % 8 == 2 ) ? 2 :
            ( i % 8 == 3 ) ? 3 :
            ( i % 8 == 4 ) ? 7 :
            ( i % 8 == 5 ) ? 15 :
            ( i % 8 == 6 ) ? 31 : 63
        );
        cur += step;
        v.push_back( cur );
    }

    return v;
}

template <typename UInt, typename SInt>
std::vector<UInt> make_zigzag_like_values_( std::size_t n )
{
    static_assert( sizeof(UInt) == sizeof(SInt), "signed/unsigned size mismatch" );

    std::vector<UInt> v;
    v.reserve( n );

    for( std::size_t i = 0; i < n; ++i ) {
        SInt const s = static_cast<SInt>(
            ( i % 2 == 0 )
                ? static_cast<SInt>( i * 3 )
                : -static_cast<SInt>( i * 5 + 1 )
        );
        v.push_back( static_cast<UInt>( s ));
    }

    return v;
}

template <typename UInt>
void expect_equal_roundtrip_(
    std::vector<UInt> const& input,
    std::vector<UInt> const& output
) {
    ASSERT_EQ( input.size(), output.size() );
    EXPECT_EQ( input, output );
}

// =================================================================================================
//     32-bit helpers
// =================================================================================================

void roundtrip_plain_32_( std::vector<std::uint32_t> const& input )
{
    auto compressed = std::vector<std::uint8_t>( pfor_bound_32( input.size() ));
    auto decoded    = std::vector<std::uint32_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_32(
        working_in.data(), working_in.size(), compressed.data()
    );

    if( input.empty() ) {
        EXPECT_EQ( used, 0u );
    } else {
        ASSERT_GT( used, 0u );
        ASSERT_LE( used, compressed.size() );
    }

    std::size_t const consumed = pfor_decode_32(
        compressed.data(), input.size(), decoded.data()
    );

    if( input.empty() ) {
        EXPECT_EQ( consumed, 0u );
    } else {
        ASSERT_GT( consumed, 0u );
        ASSERT_LE( consumed, used );
    }

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_delta_32_( std::vector<std::uint32_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_32( input.size() ));
    auto decoded    = std::vector<std::uint32_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_delta_32(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_delta_32(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_delta1_32_( std::vector<std::uint32_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_32( input.size() ));
    auto decoded    = std::vector<std::uint32_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_delta1_32(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_delta1_32(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_zigzag_32_( std::vector<std::uint32_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_32( input.size() ));
    auto decoded    = std::vector<std::uint32_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_zigzag_32(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_zigzag_32(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

// =================================================================================================
//     64-bit helpers
// =================================================================================================

void roundtrip_plain_64_( std::vector<std::uint64_t> const& input )
{
    auto compressed = std::vector<std::uint8_t>( pfor_bound_64( input.size() ));
    auto decoded    = std::vector<std::uint64_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_64(
        working_in.data(), working_in.size(), compressed.data()
    );

    if( input.empty() ) {
        EXPECT_EQ( used, 0u );
    } else {
        ASSERT_GT( used, 0u );
        ASSERT_LE( used, compressed.size() );
    }

    std::size_t const consumed = pfor_decode_64(
        compressed.data(), input.size(), decoded.data()
    );

    if( input.empty() ) {
        EXPECT_EQ( consumed, 0u );
    } else {
        ASSERT_GT( consumed, 0u );
        ASSERT_LE( consumed, used );
    }

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_delta_64_( std::vector<std::uint64_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_64( input.size() ));
    auto decoded    = std::vector<std::uint64_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_delta_64(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_delta_64(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_delta1_64_( std::vector<std::uint64_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_64( input.size() ));
    auto decoded    = std::vector<std::uint64_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_delta1_64(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_delta1_64(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

void roundtrip_zigzag_64_( std::vector<std::uint64_t> const& input )
{
    ASSERT_FALSE( input.empty() );

    auto compressed = std::vector<std::uint8_t>( pfor_bound_64( input.size() ));
    auto decoded    = std::vector<std::uint64_t>( input.size(), 0 );

    auto working_in = input;
    std::size_t const used = pfor_encode_zigzag_64(
        working_in.data(), working_in.size(), compressed.data()
    );
    ASSERT_GT( used, 0u );
    ASSERT_LE( used, compressed.size() );

    std::size_t const consumed = pfor_decode_zigzag_64(
        compressed.data(), input.size(), decoded.data()
    );
    ASSERT_GT( consumed, 0u );
    ASSERT_LE( consumed, used );

    expect_equal_roundtrip_( input, decoded );
}

// =================================================================================================
//     Tests
// =================================================================================================

TEST( Pfor, Roundtrip32BitAllFlavours )
{
    std::vector<std::vector<std::uint32_t>> plain_inputs{
        {},
        { 0u },
        { 1u, 1u, 1u, 1u, 1u },
        make_patterned_values_<std::uint32_t>( 17 ),
        make_patterned_values_<std::uint32_t>( 128 ),
        make_patterned_values_<std::uint32_t>( 257 )
    };

    std::vector<std::vector<std::uint32_t>> monotone_inputs{
        { 0u },
        { 1u, 2u, 3u, 4u, 5u },
        make_monotone_values_<std::uint32_t>( 17 ),
        make_monotone_values_<std::uint32_t>( 128 ),
        make_monotone_values_<std::uint32_t>( 257 )
    };

    std::vector<std::vector<std::uint32_t>> zigzag_inputs{
        { 0u },
        make_zigzag_like_values_<std::uint32_t, std::int32_t>( 17 ),
        make_zigzag_like_values_<std::uint32_t, std::int32_t>( 128 ),
        make_zigzag_like_values_<std::uint32_t, std::int32_t>( 257 )
    };

    for( auto const& input : plain_inputs ) {
        roundtrip_plain_32_( input );
    }
    for( auto const& input : monotone_inputs ) {
        roundtrip_delta_32_( input );
        roundtrip_delta1_32_( input );
    }
    for( auto const& input : zigzag_inputs ) {
        roundtrip_zigzag_32_( input );
    }
}

TEST( Pfor, Roundtrip64BitAllFlavours )
{
    std::vector<std::vector<std::uint64_t>> plain_inputs{
        {},
        { 0ull },
        { 1ull, 1ull, 1ull, 1ull, 1ull },
        make_patterned_values_<std::uint64_t>( 17 ),
        make_patterned_values_<std::uint64_t>( 128 ),
        make_patterned_values_<std::uint64_t>( 257 )
    };

    std::vector<std::vector<std::uint64_t>> monotone_inputs{
        { 0ull },
        { 1ull, 2ull, 3ull, 4ull, 5ull },
        make_monotone_values_<std::uint64_t>( 17 ),
        make_monotone_values_<std::uint64_t>( 128 ),
        make_monotone_values_<std::uint64_t>( 257 )
    };

    std::vector<std::vector<std::uint64_t>> zigzag_inputs{
        { 0ull },
        make_zigzag_like_values_<std::uint64_t, std::int64_t>( 17 ),
        make_zigzag_like_values_<std::uint64_t, std::int64_t>( 128 ),
        make_zigzag_like_values_<std::uint64_t, std::int64_t>( 257 )
    };

    for( auto const& input : plain_inputs ) {
        roundtrip_plain_64_( input );
    }
    for( auto const& input : monotone_inputs ) {
        roundtrip_delta_64_( input );
        roundtrip_delta1_64_( input );
    }
    for( auto const& input : zigzag_inputs ) {
        roundtrip_zigzag_64_( input );
    }
}

TEST( Pfor, ZeroLengthInputReturnsZero )
{
    std::vector<std::uint32_t> in32;
    std::vector<std::uint64_t> in64;
    std::vector<std::uint8_t>  out;

    std::vector<std::uint32_t> dec32;
    std::vector<std::uint64_t> dec64;

    EXPECT_EQ( pfor_encode_32( in32.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_32( out.data(), 0, dec32.data() ), 0u );

    EXPECT_EQ( pfor_encode_delta_32( in32.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_delta_32( out.data(), 0, dec32.data() ), 0u );

    EXPECT_EQ( pfor_encode_delta1_32( in32.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_delta1_32( out.data(), 0, dec32.data() ), 0u );

    EXPECT_EQ( pfor_encode_zigzag_32( in32.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_zigzag_32( out.data(), 0, dec32.data() ), 0u );

    EXPECT_EQ( pfor_encode_64( in64.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_64( out.data(), 0, dec64.data() ), 0u );

    EXPECT_EQ( pfor_encode_delta_64( in64.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_delta_64( out.data(), 0, dec64.data() ), 0u );

    EXPECT_EQ( pfor_encode_delta1_64( in64.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_delta1_64( out.data(), 0, dec64.data() ), 0u );

    EXPECT_EQ( pfor_encode_zigzag_64( in64.data(), 0, out.data() ), 0u );
    EXPECT_EQ( pfor_decode_zigzag_64( out.data(), 0, dec64.data() ), 0u );
}

// =================================================================================================
//     Performance Tests
// =================================================================================================

// Strictly monotonic, steps cycling 1..8; values stay well within uint32 range.
template<typename UInt>
static std::vector<UInt> make_perf_small_monotonic_( std::size_t n )
{
    std::vector<UInt> v;
    v.reserve(n);
    UInt cur = 0;
    for( std::size_t i = 0; i < n; ++i ) {
        cur += static_cast<UInt>(1 + (i % 8));
        v.push_back(cur);
    }
    return v;
}

// Strictly monotonic with large steps; deltas are multiples of 2^32 so they span the full
// 64-bit range and require all 64 bits to represent each delta value.
static std::vector<std::uint64_t> make_perf_large_monotonic_( std::size_t n )
{
    std::vector<std::uint64_t> v;
    v.reserve(n);
    std::uint64_t cur = 0;
    for( std::size_t i = 0; i < n; ++i ) {
        cur += (std::uint64_t(1) + (i % 8)) * std::uint64_t(0x1'0000'0000ULL);
        v.push_back(cur);
    }
    return v;
}

template<typename UInt, typename Encode, typename Decode, typename Bound>
static void run_perf_delta_(
    char const*               label,
    std::vector<UInt> const&  input,
    Encode                    encode,
    Decode                    decode,
    Bound                     bound,
    double                    target_s
) {
    std::size_t const n = input.size();
    auto working    = input;  // encode modifies buffer in-place, so we copy before each call
    auto compressed = std::vector<std::uint8_t>(bound(n));
    auto decoded    = std::vector<UInt>(n);

    // Encode phase: copy input back each iteration to undo in-place delta transform.
    std::size_t compressed_bytes = 0;
    std::size_t enc_iters = 0;
    double enc_elapsed = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    do {
        std::copy(input.begin(), input.end(), working.begin());
        compressed_bytes = encode(working.data(), n, compressed.data());
        ++enc_iters;
        enc_elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    } while (enc_elapsed < target_s);
    double enc_gints_s = static_cast<double>(enc_iters) * n / enc_elapsed / 1e9;
    double ratio = static_cast<double>(n * sizeof(UInt)) /
                   static_cast<double>(compressed_bytes);

    // Decode phase: compressed buffer is reused unchanged across iterations.
    std::size_t dec_iters = 0;
    double dec_elapsed = 0.0;
    t0 = std::chrono::steady_clock::now();
    do {
        decode(compressed.data(), n, decoded.data());
        ++dec_iters;
        dec_elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    } while (dec_elapsed < target_s);
    double dec_gints_s = static_cast<double>(dec_iters) * n / dec_elapsed / 1e9;

    std::printf(
        "  %-54s  enc %6.3f Gint/s  dec %6.3f Gint/s  ratio %5.2f:1\n",
        label, enc_gints_s, dec_gints_s, ratio
    );
    std::fflush(stdout);
}

// Perf test: strictly monotonic lists, delta encoding.
//
// Three variants are measured:
//   A  uint32_t values with small steps (1..8) using the 32-bit pfor implementation.
//   B  The same numeric values widened to uint64_t, using the 64-bit pfor implementation.
//      Directly comparable to A: isolates 32-bit vs 64-bit implementation overhead on
//      identical data ranges.
//   C  uint64_t values with large steps (multiples of 2^32), using the 64-bit pfor
//      implementation.  Deltas require all 64 bits, so compression ratio drops and the
//      codec must work harder per value than in B.
//
// The encode phase includes the cost of restoring the input before each call (delta encoding
// transforms the buffer in-place).  The decode phase reuses the same compressed buffer.
TEST( Pfor, DISABLED_PerformanceComparison )
{
    constexpr std::size_t N = 16'000'000;
    constexpr double      T = 10.0;

    auto const v32_small = make_perf_small_monotonic_<std::uint32_t>( N );
    // Numerically identical to v32_small but stored as uint64_t.
    auto const v64_small = std::vector<std::uint64_t>( v32_small.begin(), v32_small.end() );
    auto const v64_large = make_perf_large_monotonic_( N );

    std::printf(
        "\n=== PFor delta encode/decode performance  (N = %zu, target %.0f s per phase) ===\n\n",
        N, T
    );

    // A: uint32 baseline
    run_perf_delta_(
        "uint32  small values  pfor_delta_32",
        v32_small,
        pfor_encode_delta_32, pfor_decode_delta_32, pfor_bound_32,
        T
    );

    // B: same values, 64-bit pfor
    run_perf_delta_(
        "uint64  small values (same as uint32 range)  pfor_delta_64",
        v64_small,
        pfor_encode_delta_64, pfor_decode_delta_64, pfor_bound_64,
        T
    );

    // C: large 64-bit values, 64-bit pfor
    run_perf_delta_(
        "uint64  large values (full 64-bit range)  pfor_delta_64",
        v64_large,
        pfor_encode_delta_64, pfor_decode_delta_64, pfor_bound_64,
        T
    );

    std::printf("\n");
}

// =================================================================================================
//     Performance: all PFor wrappers
// =================================================================================================

static std::vector<std::uint32_t> make_perf_random_32_( std::size_t n )
{
    std::mt19937 rng(42);
    std::uniform_int_distribution<std::uint32_t> dist;
    std::vector<std::uint32_t> v(n);
    std::generate(v.begin(), v.end(), [&]{ return dist(rng); });
    return v;
}

static std::vector<std::uint64_t> make_perf_random_64_( std::size_t n )
{
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<std::uint64_t> dist;
    std::vector<std::uint64_t> v(n);
    std::generate(v.begin(), v.end(), [&]{ return dist(rng); });
    return v;
}

// Oscillating walk with small signed steps, stored as uint (bit-reinterpretation).
// Consecutive differences stay in [-64, 63], which zigzag-maps to [0, 127], so PFor
// achieves near-maximum compression for this variant.
static std::vector<std::uint32_t> make_perf_zigzag_32_( std::size_t n )
{
    std::vector<std::uint32_t> v(n);
    std::int32_t cur = 0;
    for( std::size_t i = 0; i < n; ++i ) {
        std::int32_t const step = static_cast<std::int32_t>( i % 128 ) - 64;
        cur += step;
        v[i] = static_cast<std::uint32_t>( cur );
    }
    return v;
}

static std::vector<std::uint64_t> make_perf_zigzag_64_( std::size_t n )
{
    std::vector<std::uint64_t> v(n);
    std::int64_t cur = 0;
    for( std::size_t i = 0; i < n; ++i ) {
        std::int64_t const step = static_cast<std::int64_t>( i % 128 ) - 64;
        cur += step;
        v[i] = static_cast<std::uint64_t>( cur );
    }
    return v;
}

// Perf test: all PFor wrappers across all four data shapes.
//
// Four data shapes are tested for each of 32-bit and 64-bit:
//   plain   - fully random integers; worst-case compression, exercises the codec with
//             incompressible bit patterns.
//   delta   - random integers sorted ascending (non-decreasing); deltas are small on
//             average (≈ max_val/N), compresses well.
//   delta1  - same sorted sequence with duplicates removed (strictly increasing); required
//             by the delta-1 codec, typical of postings lists.
//   zigzag  - oscillating walk with signed steps in [-64, 63]; the codec applies
//             delta-then-zigzag internally, so these small signed consecutive differences
//             compress extremely well.
TEST( Pfor, DISABLED_PerformancePforWrappers )
{
    constexpr std::size_t N = 16'000'000;
    constexpr double      T = 5.0;

    // --- build 32-bit datasets ---
    auto rand32   = make_perf_random_32_(N);
    auto sorted32 = rand32;
    std::sort( sorted32.begin(), sorted32.end() );
    auto unique32 = sorted32;
    unique32.erase( std::unique( unique32.begin(), unique32.end() ), unique32.end() );
    auto zigzag32 = make_perf_zigzag_32_(N);

    // --- build 64-bit datasets ---
    auto rand64   = make_perf_random_64_(N);
    auto sorted64 = rand64;
    std::sort( sorted64.begin(), sorted64.end() );
    auto unique64 = sorted64;
    unique64.erase( std::unique( unique64.begin(), unique64.end() ), unique64.end() );
    auto zigzag64 = make_perf_zigzag_64_(N);

    std::printf(
        "\n=== PFor wrapper encode/decode performance  (N = %zu, target %.0f s per phase) ===\n\n",
        N, T
    );

    std::printf("--- 32-bit ---\n");
    run_perf_delta_( "plain   random                   pfor_32",
        rand32,   pfor_encode_32,        pfor_decode_32,        pfor_bound_32, T );
    run_perf_delta_( "delta   sorted random            pfor_delta_32",
        sorted32, pfor_encode_delta_32,  pfor_decode_delta_32,  pfor_bound_32, T );
    run_perf_delta_( "delta1  sorted unique random     pfor_delta1_32",
        unique32, pfor_encode_delta1_32, pfor_decode_delta1_32, pfor_bound_32, T );
    run_perf_delta_( "zigzag  oscillating signed steps  pfor_zigzag_32",
        zigzag32, pfor_encode_zigzag_32, pfor_decode_zigzag_32, pfor_bound_32, T );

    std::printf("\n--- 64-bit ---\n");
    run_perf_delta_( "plain   random                   pfor_64",
        rand64,   pfor_encode_64,        pfor_decode_64,        pfor_bound_64, T );
    run_perf_delta_( "delta   sorted random            pfor_delta_64",
        sorted64, pfor_encode_delta_64,  pfor_decode_delta_64,  pfor_bound_64, T );
    run_perf_delta_( "delta1  sorted unique random     pfor_delta1_64",
        unique64, pfor_encode_delta1_64, pfor_decode_delta1_64, pfor_bound_64, T );
    run_perf_delta_( "zigzag  oscillating signed steps  pfor_zigzag_64",
        zigzag64, pfor_encode_zigzag_64, pfor_decode_zigzag_64, pfor_bound_64, T );

    std::printf("\n");
}
