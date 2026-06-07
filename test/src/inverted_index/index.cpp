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
 * @brief Tests for InvertedIndex (reader / accessor).
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "spear/inverted_index/builder.hpp"
#include "spear/inverted_index/format.hpp"
#include "spear/inverted_index/index.hpp"
#include "genesis/util/threading/thread_pool.hpp"

using namespace spear::inverted_index;

using Builder = InvertedIndexBuilder<>;
using Index   = InvertedIndex<>;

// =================================================================================================
//     Shared Fixture
// =================================================================================================

// Build a 5-term index with max_postings_per_term = 4 and write it to `path`.
//
// Term layout:
//   0 : empty   (no postings)
//   1 : normal  {10, 20, 30}
//   2 : normal  {100, 200}
//   3 : capped  (5 unique postings > cap of 4)
//   4 : normal  {42}
//
// TempFile unlinks the placeholder before this call, so no overwrite check is needed.
static void build_test_index( std::string const& path )
{
    Builder::Config cfg;
    cfg.max_postings_per_term = 4;
    Builder builder( 5, cfg );

    builder.add( 1, 10 );
    builder.add( 1, 20 );
    builder.add( 1, 30 );

    builder.add( 2, 100 );
    builder.add( 2, 200 );

    // Five unique values → exceeds the cap of 4; written as capped sentinel
    builder.add( 3, 1 );
    builder.add( 3, 2 );
    builder.add( 3, 3 );
    builder.add( 3, 4 );
    builder.add( 3, 5 );

    builder.add( 4, 42 );

    builder.write( path );
}

// =================================================================================================
//     Tests
// =================================================================================================

TEST( InvertedIndex, Metadata )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx;
    EXPECT_FALSE( idx.is_open() );

    idx.open( tmp.path );
    EXPECT_TRUE( idx.is_open() );
    EXPECT_EQ( idx.num_terms(), 5u );
    EXPECT_EQ( idx.max_postings_per_term(), 4u );
    EXPECT_EQ( idx.capped_sentinel(), 5u );
}

TEST( InvertedIndex, RoundTripLoadAll )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );

    // Term 0: empty
    EXPECT_EQ( idx.posting_count( 0 ), 0u );
    EXPECT_FALSE( idx.is_capped( 0 ) );
    {
        auto const [posts, status] = idx.postings( 0 );
        EXPECT_EQ( status, Index::PostingsStatus::kEmpty );
        EXPECT_TRUE( posts.empty() );
    }

    // Term 1: normal
    EXPECT_EQ( idx.posting_count( 1 ), 3u );
    EXPECT_FALSE( idx.is_capped( 1 ) );
    {
        auto const [posts, status] = idx.postings( 1 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 10, 20, 30 } ));
    }

    // Term 2: normal
    EXPECT_EQ( idx.posting_count( 2 ), 2u );
    EXPECT_FALSE( idx.is_capped( 2 ) );
    {
        auto const [posts, status] = idx.postings( 2 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 100, 200 } ));
    }

    // Term 3: capped
    EXPECT_EQ( idx.posting_count( 3 ), idx.capped_sentinel() );
    EXPECT_TRUE( idx.is_capped( 3 ) );
    {
        auto const [posts, status] = idx.postings( 3 );
        EXPECT_EQ( status, Index::PostingsStatus::kCapped );
        EXPECT_TRUE( posts.empty() );
    }

    // Term 4: normal
    EXPECT_EQ( idx.posting_count( 4 ), 1u );
    EXPECT_FALSE( idx.is_capped( 4 ) );
    {
        auto const [posts, status] = idx.postings( 4 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 42 } ));
    }
}

TEST( InvertedIndex, RoundTripPread )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kPread );

    // Term 0: empty
    EXPECT_EQ( idx.posting_count( 0 ), 0u );
    EXPECT_FALSE( idx.is_capped( 0 ) );
    {
        auto const [posts, status] = idx.postings( 0 );
        EXPECT_EQ( status, Index::PostingsStatus::kEmpty );
        EXPECT_TRUE( posts.empty() );
    }

    // Term 1: normal
    EXPECT_EQ( idx.posting_count( 1 ), 3u );
    EXPECT_FALSE( idx.is_capped( 1 ) );
    {
        auto const [posts, status] = idx.postings( 1 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 10, 20, 30 } ));
    }

    // Term 2: normal
    EXPECT_EQ( idx.posting_count( 2 ), 2u );
    EXPECT_FALSE( idx.is_capped( 2 ) );
    {
        auto const [posts, status] = idx.postings( 2 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 100, 200 } ));
    }

    // Term 3: capped
    EXPECT_EQ( idx.posting_count( 3 ), idx.capped_sentinel() );
    EXPECT_TRUE( idx.is_capped( 3 ) );
    {
        auto const [posts, status] = idx.postings( 3 );
        EXPECT_EQ( status, Index::PostingsStatus::kCapped );
        EXPECT_TRUE( posts.empty() );
    }

    // Term 4: normal
    EXPECT_EQ( idx.posting_count( 4 ), 1u );
    EXPECT_FALSE( idx.is_capped( 4 ) );
    {
        auto const [posts, status] = idx.postings( 4 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        EXPECT_EQ( posts, ( std::vector<std::uint64_t>{ 42 } ));
    }
}

TEST( InvertedIndex, OutOfRange )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx;
    idx.open( tmp.path );

    EXPECT_THROW( (void) idx.postings( 5 ),       std::out_of_range );
    EXPECT_THROW( (void) idx.posting_count( 5 ),  std::out_of_range );
    EXPECT_THROW( (void) idx.is_capped( 5 ),      std::out_of_range );
}

TEST( InvertedIndex, BadMagic )
{
    TempFile tmp;
    build_test_index( tmp.path );

    // Corrupt the first byte of the footer (start of the magic field)
    {
        std::fstream f( tmp.path, std::ios::in | std::ios::out | std::ios::binary );
        ASSERT_TRUE( f.is_open() );
        f.seekp( -static_cast<std::streamoff>( sizeof( InvertedIndexFooter )), std::ios::end );
        char const bad = '\xFF';
        f.write( &bad, 1 );
    }

    Index idx;
    EXPECT_THROW( idx.open( tmp.path ), std::runtime_error );
}

TEST( InvertedIndex, WrongPositionBits )
{
    TempFile tmp;
    build_test_index( tmp.path );  // written with uint64_t (64-bit positions)

    InvertedIndex<std::uint32_t> idx;
    EXPECT_THROW( idx.open( tmp.path ), std::runtime_error );
}

TEST( InvertedIndex, MoveConstructor )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx1;
    idx1.open( tmp.path, Index::OpenMode::kPread );
    EXPECT_TRUE( idx1.is_open() );

    Index idx2 = std::move( idx1 );

    EXPECT_FALSE( idx1.is_open() );
    EXPECT_TRUE( idx2.is_open() );
    EXPECT_EQ( idx2.postings( 1 ).first, ( std::vector<std::uint64_t>{ 10, 20, 30 } ));
}

TEST( InvertedIndex, MoveAssignment )
{
    TempFile tmp;
    build_test_index( tmp.path );

    Index idx1;
    idx1.open( tmp.path, Index::OpenMode::kPread );

    Index idx2;
    idx2 = std::move( idx1 );

    EXPECT_FALSE( idx1.is_open() );
    EXPECT_TRUE( idx2.is_open() );
    EXPECT_EQ( idx2.postings( 4 ).first, ( std::vector<std::uint64_t>{ 42 } ));
}

TEST( InvertedIndex, Reopen )
{
    // Two independent temp files, two independent builders
    TempFile tmp1;
    TempFile tmp2;

    {
        Builder b1( 3 );
        b1.add( 0, 5 );
        b1.add( 2, 77 );
        b1.write( tmp1.path );
    }
    {
        Builder b2( 2 );
        b2.add( 0, 99 );
        b2.write( tmp2.path );
    }

    Index idx;

    idx.open( tmp1.path, Index::OpenMode::kPread );
    EXPECT_EQ( idx.num_terms(), 3u );
    EXPECT_EQ( idx.postings( 0 ).first, ( std::vector<std::uint64_t>{ 5 } ));

    // Reopen with a different file — old fd must be released cleanly
    idx.open( tmp2.path, Index::OpenMode::kPread );
    EXPECT_EQ( idx.num_terms(), 2u );
    EXPECT_EQ( idx.postings( 0 ).first, ( std::vector<std::uint64_t>{ 99 } ));
}

// =================================================================================================
//     No-Cap Tests
// =================================================================================================

// Build a no-cap index: 3 terms, term 0 has 2000 postings (above the default cap of 1024),
// term 1 has 5 postings, term 2 is empty. Verify round-trip in both open modes.
static void build_nocap_index( std::string const& path )
{
    Builder::Config cfg;
    cfg.max_postings_per_term = 0;   // no capping
    Builder builder( 3, cfg );

    for( std::uint64_t p = 0; p < 2000; ++p ) {
        builder.add( 0, p );
    }
    builder.add( 1, 10 );
    builder.add( 1, 20 );
    builder.add( 1, 30 );
    builder.add( 1, 40 );
    builder.add( 1, 50 );
    // term 2: empty

    builder.write( path );
}

TEST( InvertedIndex, NoCapRoundTripLoadAll )
{
    TempFile tmp;
    build_nocap_index( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );

    // Term 0: 2000 postings — would be capped under the default limit of 1024
    EXPECT_EQ( idx.posting_count( 0 ), 2000u );
    EXPECT_FALSE( idx.is_capped( 0 ) );
    {
        auto const [posts, status] = idx.postings( 0 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        ASSERT_EQ( posts.size(), 2000u );
        EXPECT_EQ( posts.front(), 0u );
        EXPECT_EQ( posts.back(), 1999u );
    }

    // Term 1: normal
    EXPECT_EQ( idx.posting_count( 1 ), 5u );
    EXPECT_FALSE( idx.is_capped( 1 ) );
    EXPECT_EQ( idx.postings( 1 ).first, ( std::vector<std::uint64_t>{ 10, 20, 30, 40, 50 } ));

    // Term 2: empty
    EXPECT_EQ( idx.posting_count( 2 ), 0u );
    EXPECT_FALSE( idx.is_capped( 2 ) );
    EXPECT_EQ( idx.postings( 2 ).second, Index::PostingsStatus::kEmpty );
}

TEST( InvertedIndex, NoCapRoundTripPread )
{
    TempFile tmp;
    build_nocap_index( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kPread );

    // Term 0: 2000 postings — would be capped under the default limit of 1024
    EXPECT_EQ( idx.posting_count( 0 ), 2000u );
    EXPECT_FALSE( idx.is_capped( 0 ) );
    {
        auto const [posts, status] = idx.postings( 0 );
        EXPECT_EQ( status, Index::PostingsStatus::kFound );
        ASSERT_EQ( posts.size(), 2000u );
        EXPECT_EQ( posts.front(), 0u );
        EXPECT_EQ( posts.back(), 1999u );
    }

    // Term 1: normal
    EXPECT_EQ( idx.posting_count( 1 ), 5u );
    EXPECT_FALSE( idx.is_capped( 1 ) );
    EXPECT_EQ( idx.postings( 1 ).first, ( std::vector<std::uint64_t>{ 10, 20, 30, 40, 50 } ));

    // Term 2: empty
    EXPECT_EQ( idx.posting_count( 2 ), 0u );
    EXPECT_FALSE( idx.is_capped( 2 ) );
    EXPECT_EQ( idx.postings( 2 ).second, Index::PostingsStatus::kEmpty );
}

// =================================================================================================
//     Performance Test
// =================================================================================================

// ~1 GB of uncompressed uint64_t postings (128M values × 8 bytes), distributed randomly across
// 1M terms. Positions are added in strictly increasing global order so delta-1 encoding is valid.
// No capping (max_postings_per_term = 65534, the largest safe value for uint16_t StoredCountT).
//
// Metrics printed (no hard thresholds — they vary by hardware):
//   build time, write time, file size + compression ratio,
//   kLoadAll open + full-scan query time, kPread full-scan query time.
//
// File is written to $HOME/spear_perf_test.idx and NOT deleted so the result can be inspected.
TEST( InvertedIndex, DISABLED_LargeScalePerf )
{
    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::duration<double, std::milli>;

    constexpr std::size_t   num_terms    = 1'000'000;
    constexpr std::uint64_t num_postings = 128'000'000;   // × 8 bytes = 1 GB uncompressed
    constexpr double        uncompressed_gb = static_cast<double>( num_postings ) * 8.0 / 1e9;

    char const* home = std::getenv( "HOME" );
    std::string const path = std::string( home ? home : "/tmp" ) + "/spear_perf_test.idx";
    std::cout << "\n[PerfTest] output: " << path << "\n";

    // ---- Build ----
    Builder::Config cfg;
    cfg.max_postings_per_term = 0;   // no capping
    Builder builder( num_terms, cfg );

    std::mt19937_64 rng( 42 );
    std::uniform_int_distribution<std::size_t> term_dist( 0, num_terms - 1 );

    auto t0 = Clock::now();
    for( std::uint64_t p = 0; p < num_postings; ++p ) {
        builder.add( term_dist( rng ), p );
    }
    double const build_ms = Ms( Clock::now() - t0 ).count();
    std::cout << "[PerfTest] build:  " << build_ms / 1e3 << " s"
              << "  (" << num_postings / ( build_ms / 1e3 ) / 1e6 << " M postings/s)\n";

    // ---- Write ----
    std::remove( path.c_str() );  // delete previous run's output if present
    auto t1 = Clock::now();
    builder.write( path );
    double const write_ms = Ms( Clock::now() - t1 ).count();

    // File size
    std::FILE* fsize = std::fopen( path.c_str(), "rb" );
    ASSERT_TRUE( fsize != nullptr );
    std::fseek( fsize, 0, SEEK_END );
    double const file_bytes = static_cast<double>( std::ftell( fsize ) );
    std::fclose( fsize );
    double const file_gb = file_bytes / 1e9;

    std::cout << "[PerfTest] write:  " << write_ms / 1e3 << " s"
              << "  (" << file_gb / 1e9 * 1e9 / 1e6 << " GB written)\n";
    std::cout << "[PerfTest] size:   " << file_gb << " GB"
              << "  (ratio " << uncompressed_gb / file_gb << "x vs "
              << uncompressed_gb << " GB uncompressed)\n";

    // ---- kLoadAll: open + full scan ----
    {
        Index idx;
        auto t2 = Clock::now();
        idx.open( path, Index::OpenMode::kLoadAll );
        double const open_ms = Ms( Clock::now() - t2 ).count();

        std::vector<std::uint64_t> buf;
        std::uint64_t checksum = 0;
        std::uint64_t total_postings = 0;
        auto t3 = Clock::now();
        for( std::size_t t = 0; t < num_terms; ++t ) {
            idx.postings( t, buf );
            total_postings += buf.size();
            if( !buf.empty() ) {
                checksum ^= buf.front() ^ buf.back();
            }
        }
        double const scan_ms = Ms( Clock::now() - t3 ).count();

        std::cout << "[PerfTest] kLoadAll open:  " << open_ms / 1e3 << " s\n";
        std::cout << "[PerfTest] kLoadAll scan:  " << scan_ms / 1e3 << " s"
                  << "  (" << num_terms / ( scan_ms / 1e3 ) / 1e6 << " M terms/s)"
                  << "  total_postings=" << total_postings
                  << "  checksum=" << checksum << "\n";
    }

    // ---- kPread: open + full scan ----
    {
        Index idx;
        auto t4 = Clock::now();
        idx.open( path, Index::OpenMode::kPread );
        double const open_ms = Ms( Clock::now() - t4 ).count();

        std::vector<std::uint64_t> buf;
        std::uint64_t checksum = 0;
        std::uint64_t total_postings = 0;
        auto t5 = Clock::now();
        for( std::size_t t = 0; t < num_terms; ++t ) {
            idx.postings( t, buf );
            total_postings += buf.size();
            if( !buf.empty() ) {
                checksum ^= buf.front() ^ buf.back();
            }
        }
        double const scan_ms = Ms( Clock::now() - t5 ).count();

        std::cout << "[PerfTest] kPread  open:  " << open_ms / 1e3 << " s\n";
        std::cout << "[PerfTest] kPread  scan:  " << scan_ms / 1e3 << " s"
                  << "  (" << num_terms / ( scan_ms / 1e3 ) / 1e6 << " M terms/s)"
                  << "  total_postings=" << total_postings
                  << "  checksum=" << checksum << "\n";
    }
}

// =================================================================================================
//     Concurrent Stress Tests
// =================================================================================================

// P = number of positions each thread inserts per term. With T threads and P positions,
// each term receives T*P unique postings in total. Raise for a longer soak test; keep small for CI.
static constexpr std::size_t kStressP = 16;

// Thread t inserts positions t, t+T, t+2T, …, t+(P-1)*T per term (T = num_threads, stride T).
// All threads together cover exactly {0, 1, …, T*P-1} with no duplicates, so after a round-trip
// each term must have exactly T*P postings equal to {0, 1, …, T*P-1}. No reference map needed.
static void run_stress_nocap_(
    std::size_t num_threads,
    std::size_t num_terms,
    std::size_t pending_capacity,
    std::size_t P
) {
    using namespace genesis::util::threading;
    TempFile tmp;

    Builder::Config cfg;
    cfg.max_postings_per_term = 0;
    cfg.pending_capacity      = pending_capacity;
    Builder builder( num_terms, cfg );

    auto pool = std::make_shared<ThreadPool>( num_threads );
    std::vector<ProactiveFuture<bool>> futures;
    futures.reserve( num_threads );

    for( std::size_t t = 0; t < num_threads; ++t ) {
        futures.emplace_back( pool->enqueue_and_retrieve(
            [&builder, t, num_threads, num_terms, P]() -> bool
            {
                std::vector<std::size_t> order( num_terms );
                std::iota( order.begin(), order.end(), 0u );
                std::mt19937 rng( static_cast<std::uint32_t>( t ));
                for( std::size_t k = 0; k < P; ++k ) {
                    std::shuffle( order.begin(), order.end(), rng );
                    for( std::size_t term : order ) {
                        builder.add( term, static_cast<std::uint64_t>( t + k * num_threads ));
                    }
                }
                return true;
            }
        ));
    }
    for( auto& f : futures ) {
        f.get();
    }

    builder.write( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );
    ASSERT_EQ( idx.num_terms(), num_terms );

    std::vector<std::uint64_t> expected( num_threads * P );
    std::iota( expected.begin(), expected.end(), std::uint64_t{0} );

    for( std::size_t term = 0; term < num_terms; ++term ) {
        ASSERT_FALSE( idx.is_capped( term ))
            << "term=" << term << " threads=" << num_threads << " terms=" << num_terms;
        ASSERT_EQ( idx.posting_count( term ), num_threads * P )
            << "term=" << term << " threads=" << num_threads << " terms=" << num_terms;
        EXPECT_EQ( idx.postings( term ).first, expected )
            << "term=" << term << " threads=" << num_threads << " terms=" << num_terms;
    }
}

// Terms [0, num_terms/2) are "light": each thread inserts one position (= t).
// Total unique per light term = T; with cap = T they are guaranteed not to be capped (T <= T).
// Terms [num_terms/2, num_terms) are "heavy": each thread inserts t and t+T.
// Total unique per heavy term = 2T; with cap = T, 2T > T guarantees capping.
static void run_stress_capped_(
    std::size_t num_threads,
    std::size_t num_terms,
    std::size_t pending_capacity
) {
    using namespace genesis::util::threading;
    std::size_t const half = num_terms / 2;
    if( half == 0 ) {
        return;
    }

    TempFile tmp;
    Builder::Config cfg;
    cfg.max_postings_per_term = static_cast<std::uint16_t>( num_threads );
    cfg.pending_capacity      = pending_capacity;
    Builder builder( num_terms, cfg );

    auto pool = std::make_shared<ThreadPool>( num_threads );
    std::vector<ProactiveFuture<bool>> futures;
    futures.reserve( num_threads );

    for( std::size_t t = 0; t < num_threads; ++t ) {
        futures.emplace_back( pool->enqueue_and_retrieve(
            [&builder, t, num_threads, num_terms, half]() -> bool
            {
                std::vector<std::size_t> order( num_terms );
                std::iota( order.begin(), order.end(), 0u );
                std::mt19937 rng( static_cast<std::uint32_t>( t + 42 ));
                // Pass 0: insert position t into every term (light + heavy).
                std::shuffle( order.begin(), order.end(), rng );
                for( std::size_t term : order ) {
                    builder.add( term, static_cast<std::uint64_t>( t ));
                }
                // Pass 1: insert position t+T into heavy terms only, interleaved with other threads.
                std::shuffle( order.begin(), order.end(), rng );
                for( std::size_t term : order ) {
                    if( term >= half ) {
                        builder.add( term, static_cast<std::uint64_t>( t + num_threads ));
                    }
                }
                return true;
            }
        ));
    }
    for( auto& f : futures ) {
        f.get();
    }

    builder.write( tmp.path );

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );
    ASSERT_EQ( idx.num_terms(), num_terms );

    std::vector<std::uint64_t> expected_light( num_threads );
    std::iota( expected_light.begin(), expected_light.end(), std::uint64_t{0} );

    for( std::size_t term = 0; term < half; ++term ) {
        EXPECT_FALSE( idx.is_capped( term ))
            << "light term=" << term << " threads=" << num_threads;
        EXPECT_EQ( idx.posting_count( term ), num_threads )
            << "light term=" << term;
        EXPECT_EQ( idx.postings( term ).first, expected_light )
            << "light term=" << term;
    }
    for( std::size_t term = half; term < num_terms; ++term ) {
        EXPECT_TRUE( idx.is_capped( term ))
            << "heavy term=" << term << " threads=" << num_threads;
    }
}

TEST( InvertedIndex, StressNoCap )
{
    size_t const num_runs = 16;
    for( size_t i = 0; i < num_runs; ++i ) {
        for( std::size_t num_threads : { 1u, 2u, 4u, 8u } ) {
            for( std::size_t num_terms : { 1u, 64u, 1024u } ) {
                for( std::size_t pending_cap : { 2u, 16u } ) {
                    run_stress_nocap_( num_threads, num_terms, pending_cap, kStressP );
                }
            }
        }
    }
}

TEST( InvertedIndex, StressCapped )
{
    size_t const num_runs = 16;
    for( size_t i = 0; i < num_runs; ++i ) {
        for( std::size_t num_threads : { 1u, 2u, 4u, 8u } ) {
            for( std::size_t num_terms : { 2u, 64u, 1024u } ) {
                for( std::size_t pending_cap : { 2u, 16u } ) {
                    run_stress_capped_( num_threads, num_terms, pending_cap );
                }
            }
        }
    }
}

// =================================================================================================
//     Concurrent Performance Test
// =================================================================================================

// Large-scale concurrent build + write + read-back performance test.
// Disabled by default; enable with --gtest_also_run_disabled_tests.
// Raise P and/or num_terms for a longer soak; lower for a quick smoke run.
// File written to $HOME/spear_stress_perf.idx and NOT deleted for post-run inspection.
TEST( InvertedIndex, DISABLED_StressPerf )
{
    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::duration<double, std::milli>;
    using namespace genesis::util::threading;

    constexpr std::size_t   num_terms   = 1'000'000;
    constexpr std::size_t   num_threads = 8;
    constexpr std::size_t   P           = 128;   // positions per thread per term; raise for longer run

    char const* home = std::getenv( "HOME" );
    std::string const path = std::string( home ? home : "/tmp" ) + "/spear_stress_perf.idx";
    std::uint64_t const total_inserts =
        static_cast<std::uint64_t>( num_terms ) * num_threads * P;

    std::cout << "\n[StressPerfTest] output: " << path << "\n";
    std::cout << "[StressPerfTest] terms = "          << num_terms
              << "  threads = "                       << num_threads
              << "  positions per thread per term = " << P
              << "  total_inserts = "                 << total_inserts / 1'000'000 << "M\n";

    // ---- Concurrent build ----
    Builder::Config cfg;
    cfg.max_postings_per_term = 0;
    Builder builder( num_terms, cfg );

    auto pool = std::make_shared<ThreadPool>( num_threads );
    std::vector<ProactiveFuture<bool>> futures;
    futures.reserve( num_threads );

    auto t0 = Clock::now();
    for( std::size_t t = 0; t < num_threads; ++t ) {
        futures.emplace_back( pool->enqueue_and_retrieve(
            [&builder, t]() -> bool
            {
                std::vector<std::size_t> order( num_terms );
                std::iota( order.begin(), order.end(), 0u );
                std::mt19937 rng( static_cast<std::uint32_t>( t ));
                for( std::size_t k = 0; k < P; ++k ) {
                    std::shuffle( order.begin(), order.end(), rng );
                    for( std::size_t term : order ) {
                        builder.add( term, static_cast<std::uint64_t>( t + k * num_threads ));
                    }
                }
                return true;
            }
        ));
    }
    for( auto& f : futures ) {
        f.get();
    }
    double const build_ms = Ms( Clock::now() - t0 ).count();
    std::cout << "[StressPerfTest] build:  " << build_ms / 1e3 << " s"
              << "  (" << total_inserts / ( build_ms / 1e3 ) / 1e6 << " M inserts/s)\n";

    // ---- Write ----
    std::remove( path.c_str() );  // delete previous run's output if present
    auto t1 = Clock::now();
    builder.write( path );
    double const write_ms = Ms( Clock::now() - t1 ).count();
    std::FILE* fsz = std::fopen( path.c_str(), "rb" );
    ASSERT_TRUE( fsz != nullptr );
    std::fseek( fsz, 0, SEEK_END );
    double const file_gb = static_cast<double>( std::ftell( fsz )) / 1e9;
    std::fclose( fsz );
    std::cout << "[StressPerfTest] write:  " << write_ms / 1e3 << " s"
              << "  file=" << file_gb << " GB\n";

    // ---- kLoadAll: open + full scan ----
    {
        Index idx;
        auto t2 = Clock::now();
        idx.open( path, Index::OpenMode::kLoadAll );
        double const open_ms = Ms( Clock::now() - t2 ).count();

        std::vector<std::uint64_t> buf;
        std::uint64_t checksum = 0;
        std::uint64_t total_postings = 0;
        auto t3 = Clock::now();
        for( std::size_t term = 0; term < num_terms; ++term ) {
            idx.postings( term, buf );
            total_postings += buf.size();
            if( !buf.empty() ) {
                checksum ^= buf.front() ^ buf.back();
            }
        }
        double const scan_ms = Ms( Clock::now() - t3 ).count();
        std::cout << "[StressPerfTest] kLoadAll open: " << open_ms / 1e3 << " s\n";
        std::cout << "[StressPerfTest] kLoadAll scan: " << scan_ms / 1e3 << " s"
                  << "  (" << num_terms / ( scan_ms / 1e3 ) / 1e6 << " M terms/s)"
                  << "  total_postings=" << total_postings
                  << "  checksum=" << checksum << "\n";
    }

    // ---- kPread: open + full scan ----
    {
        Index idx;
        auto t4 = Clock::now();
        idx.open( path, Index::OpenMode::kPread );
        double const open_ms = Ms( Clock::now() - t4 ).count();

        std::vector<std::uint64_t> buf;
        std::uint64_t checksum = 0;
        std::uint64_t total_postings = 0;
        auto t5 = Clock::now();
        for( std::size_t term = 0; term < num_terms; ++term ) {
            idx.postings( term, buf );
            total_postings += buf.size();
            if( !buf.empty() ) {
                checksum ^= buf.front() ^ buf.back();
            }
        }
        double const scan_ms = Ms( Clock::now() - t5 ).count();
        std::cout << "[StressPerfTest] kPread  open: " << open_ms / 1e3 << " s\n";
        std::cout << "[StressPerfTest] kPread  scan: " << scan_ms / 1e3 << " s"
                  << "  (" << num_terms / ( scan_ms / 1e3 ) / 1e6 << " M terms/s)"
                  << "  total_postings=" << total_postings
                  << "  checksum=" << checksum << "\n";
    }
}
