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
 * @brief Tests for TermPostings and HitCollector.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/inverted_index/builder.hpp"
#include "spear/inverted_index/hit_collector.hpp"
#include "spear/inverted_index/index.hpp"
#include "spear/inverted_index/term_postings.hpp"

#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using namespace spear::inverted_index;

// Convenience aliases for the default (uint64_t) instantiation.
using TP       = TermPostings<std::uint64_t>;
using HC       = HitCollector<std::uint64_t, 16>;
using Interval = HC::HitInterval;

// =================================================================================================
//     Term Postings tests
// =================================================================================================

TEST( TermPostings, DirectAddVectorEmptyReturnsEmpty )
{
    TP tp( 4 );
    auto status = tp.add( std::vector<std::uint64_t>{} );
    EXPECT_EQ( status, TP::PostingsStatus::kEmpty );
    EXPECT_EQ( tp.list_count(), 0u );
    EXPECT_EQ( tp.stats().empty, 1u );
    EXPECT_EQ( tp.stats().found, 0u );
}

TEST( TermPostings, DirectAddVectorNonEmptyReturnsFound )
{
    TP tp( 4 );
    auto status = tp.add( std::vector<std::uint64_t>{ 1, 2, 3 } );
    EXPECT_EQ( status, TP::PostingsStatus::kFound );
    EXPECT_EQ( tp.list_count(), 1u );
    EXPECT_EQ( tp.stats().found, 1u );

    auto sp = tp.list_at( 0 );
    ASSERT_EQ( sp.size(), 3u );
    EXPECT_EQ( sp[0], 1u );
    EXPECT_EQ( sp[1], 2u );
    EXPECT_EQ( sp[2], 3u );
}

TEST( TermPostings, DirectAddSpanNonEmptyReturnsFound )
{
    TP tp( 4 );
    std::vector<std::uint64_t> src{ 10, 20, 30 };
    auto status = tp.add( std::span<std::uint64_t const>{ src } );
    EXPECT_EQ( status, TP::PostingsStatus::kFound );
    EXPECT_EQ( tp.list_count(), 1u );

    auto sp = tp.list_at( 0 );
    ASSERT_EQ( sp.size(), 3u );
    EXPECT_EQ( sp[0], 10u );
}

TEST( TermPostings, DirectAddExceedsCapacityThrows )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1 } );
    tp.add( std::vector<std::uint64_t>{ 2 } );
    EXPECT_THROW( tp.add( std::vector<std::uint64_t>{ 3 } ), std::runtime_error );
}

TEST( TermPostings, ClearResetsCountAndStats )
{
    TP tp( 4 );
    tp.add( std::vector<std::uint64_t>{ 1, 2 } );
    tp.add( std::vector<std::uint64_t>{} );
    EXPECT_EQ( tp.list_count(), 1u );
    EXPECT_EQ( tp.stats().found, 1u );
    EXPECT_EQ( tp.stats().empty, 1u );

    tp.clear();
    EXPECT_EQ( tp.list_count(), 0u );
    EXPECT_EQ( tp.stats().found, 0u );
    EXPECT_EQ( tp.stats().empty, 0u );

    // Verify buffers are reused (capacity not shrunk) by adding again
    tp.add( std::vector<std::uint64_t>{ 5, 6, 7 } );
    EXPECT_EQ( tp.list_count(), 1u );
}

TEST( TermPostings, StatsViaIndex )
{
    // Build a small index: term 0 empty, term 1 normal, term 2 capped.
    TempFile tmp;
    {
        InvertedIndexBuilder<>::Config cfg;
        cfg.max_postings_per_term = 3;
        InvertedIndexBuilder<> builder( 3, cfg );
        builder.add( 1, 10 );
        builder.add( 1, 20 );
        builder.add( 2, 1 );
        builder.add( 2, 2 );
        builder.add( 2, 3 );
        builder.add( 2, 4 ); // 4 values > cap of 3 → capped
        builder.write( tmp.path );
    }

    InvertedIndex<> idx;
    idx.open( tmp.path, InvertedIndex<>::OpenMode::kLoadAll );

    TP tp( 3 );
    tp.add( idx, 0 ); // kEmpty
    tp.add( idx, 1 ); // kFound
    tp.add( idx, 2 ); // kHardCapped

    EXPECT_EQ( tp.stats().empty,       1u );
    EXPECT_EQ( tp.stats().found,       1u );
    EXPECT_EQ( tp.stats().hard_capped, 1u );
    EXPECT_EQ( tp.list_count(), 1u );
}

TEST( TermPostings, SoftCap )
{
    // Build a small index with no build-time cap: term 0 has 1 posting, term 1 has 3.
    TempFile tmp;
    {
        InvertedIndexBuilder<>::Config cfg;
        InvertedIndexBuilder<> builder( 2, cfg );
        builder.add( 0, 42 );
        builder.add( 1, 10 );
        builder.add( 1, 20 );
        builder.add( 1, 30 );
        builder.write( tmp.path );
    }

    InvertedIndex<> idx;
    idx.open( tmp.path, InvertedIndex<>::OpenMode::kLoadAll );

    // Default: no soft cap, both terms found.
    {
        TP tp( 2 );
        EXPECT_EQ( tp.max_posting_list_length(), 0u );
        EXPECT_EQ( tp.add( idx, 0 ), TermPostings<std::uint64_t>::PostingsStatus::kFound );
        EXPECT_EQ( tp.add( idx, 1 ), TermPostings<std::uint64_t>::PostingsStatus::kFound );
        EXPECT_EQ( tp.stats().found,       2u );
        EXPECT_EQ( tp.stats().soft_capped, 0u );
        EXPECT_EQ( tp.list_count(), 2u );
    }

    // Soft cap of 1 via constructor: term 1 (3 postings) is soft-capped, term 0 (1) is not.
    {
        TP tp( 2, /* max_posting_list_length */ 1 );
        EXPECT_EQ( tp.max_posting_list_length(), 1u );
        EXPECT_EQ( tp.add( idx, 0 ), TermPostings<std::uint64_t>::PostingsStatus::kFound );
        EXPECT_EQ( tp.add( idx, 1 ), TermPostings<std::uint64_t>::PostingsStatus::kSoftCapped );
        EXPECT_EQ( tp.stats().found,       1u );
        EXPECT_EQ( tp.stats().soft_capped, 1u );
        EXPECT_EQ( tp.list_count(), 1u );
        auto const list0 = tp.list_at( 0 );
        EXPECT_EQ( ( std::vector<std::uint64_t>( list0.begin(), list0.end() ) ),
                   ( std::vector<std::uint64_t>{ 42 } ));
    }

    // Soft cap raised via setter: term 1 is found again.
    {
        TP tp( 2, /* max_posting_list_length */ 1 );
        tp.set_max_posting_list_length( 3 );
        EXPECT_EQ( tp.max_posting_list_length(), 3u );
        EXPECT_EQ( tp.add( idx, 0 ), TermPostings<std::uint64_t>::PostingsStatus::kFound );
        EXPECT_EQ( tp.add( idx, 1 ), TermPostings<std::uint64_t>::PostingsStatus::kFound );
        EXPECT_EQ( tp.stats().found,       2u );
        EXPECT_EQ( tp.stats().soft_capped, 0u );
        EXPECT_EQ( tp.list_count(), 2u );
    }
}

// =================================================================================================
//     HitCollector guard tests
// =================================================================================================

TEST( HitCollector, ThrowsOnMZero )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1 } );
    HC hc;
    EXPECT_THROW( hc.query( tp, 5, 0 ), std::invalid_argument );
}

TEST( HitCollector, ThrowsOnLGeRingCap )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1 } );
    HC hc; // RingCap = 16
    EXPECT_THROW( hc.query( tp, 16, 1 ), std::invalid_argument );
}

TEST( HitCollector, ThrowsOnTooManyLists )
{
    // Build a TermPostings with 257 non-empty lists to exceed the ListBitset capacity.
    TP tp( 257 );
    for( std::size_t i = 0; i < 257; ++i ) {
        tp.add( std::vector<std::uint64_t>{ i + 1 } );
    }
    HC hc;
    EXPECT_THROW( hc.query( tp, 5, 1 ), std::invalid_argument );
}

TEST( HitCollector, EmptyOnMGreaterThanListCount )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3 } );
    HC hc;
    auto result = hc.query( tp, 5, 2 ); // M=2 but only 1 list
    EXPECT_TRUE( result.empty() );
}

TEST( HitCollector, EmptyOnEmptyCollection )
{
    TP tp( 4 ); // no add() calls
    HC hc;
    auto result = hc.query( tp, 5, 1 );
    EXPECT_TRUE( result.empty() );
}

// =================================================================================================
//     Hit Collector functional tests
// =================================================================================================

// Single list, M=1: the algorithm stays in the valid region from first to last value.
// Gaps larger than L do not break the region because the single list is always in the ring.
TEST( HitCollector, SingleListM1 )
{
    TP tp( 1 );
    tp.add( std::vector<std::uint64_t>{ 10, 20, 30 } );
    HC hc;
    auto r = hc.query( tp, 5, 1 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  10u );
    EXPECT_EQ( r[0].right, 30u );
}

// Two lists that never share a window → no intervals.
TEST( HitCollector, TwoListsNoSharedWindow )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3 } );
    tp.add( std::vector<std::uint64_t>{ 20, 21, 22 } ); // gap 17 >> L=5
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    EXPECT_TRUE( r.empty() );
}

// Two lists sharing a window → one interval spanning from first co-occurring front to last back.
TEST( HitCollector, TwoListsOneInterval )
{
    // list0={1,5}, list1={3,7}, L=5, M=2
    // At v=3: both in ring [1..3] → valid, start=1
    // At v=7: evict 1 (7-1=6>5), ring=[3,5,7], both lists → valid, end=7
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 5 } );
    tp.add( std::vector<std::uint64_t>{ 3, 7 } );
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 7u );
}

// Dense overlap collapses into one maximal interval.
TEST( HitCollector, MaximalIntervalDenseOverlap )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3, 4, 5, 6, 7 } );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3, 4, 5, 6, 7 } );
    HC hc;
    auto r = hc.query( tp, 3, 2 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 7u );
}

// Wide gap (>> L) between two shared clusters produces two intervals.
// list0={1,20}, list1={3,18}, L=5, M=2:
//   cluster A: v=3 brings both lists into ring → [1,3]; v=18 evicts both → exit
//   cluster B: v=20 brings both lists back → [18,20]
TEST( HitCollector, TwoIntervalsWideGap )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 20 } );
    tp.add( std::vector<std::uint64_t>{ 3, 18 } );
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    ASSERT_EQ( r.size(), 2u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 3u );
    EXPECT_EQ( r[1].left,  18u );
    EXPECT_EQ( r[1].right, 20u );
}

// Narrow gap (< L) where the hit count drops below M and then recovers → two intervals.
// list0={1,10}, list1={4,7}, L=5, M=2:
//   v=4: both in ring [1,4] → valid [1,4]
//   v=7: evict 1 (7-1=6>5), ring=[4,7], only list1 → count=1 < 2 → emit [1,4]
//   v=10: ring=[7,10], both lists → valid [7,10] → emit [7,10]
TEST( HitCollector, TwoIntervalsNarrowGapCountDrop )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 10 } );
    tp.add( std::vector<std::uint64_t>{ 4,  7 } );
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    ASSERT_EQ( r.size(), 2u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 4u );
    EXPECT_EQ( r[1].left,  7u );
    EXPECT_EQ( r[1].right, 10u );
}

// Same value arriving from multiple lists must be deduplicated via the tail-check path
// and counted as distinct lists, not as a single list repeated.
TEST( HitCollector, DedupSameValueMultipleLists )
{
    TP tp( 3 );
    tp.add( std::vector<std::uint64_t>{ 5 } );
    tp.add( std::vector<std::uint64_t>{ 5 } );
    tp.add( std::vector<std::uint64_t>{ 5 } );
    HC hc;
    auto r = hc.query( tp, 5, 3 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  5u );
    EXPECT_EQ( r[0].right, 5u );
}

// Values exactly L apart are within the same window (val - front == L, condition is <=).
TEST( HitCollector, WindowBoundaryExactlyL )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1 } );
    tp.add( std::vector<std::uint64_t>{ 6 } ); // 6-1=5 == L=5 → same window
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 6u );
}

// Values L+1 apart are just outside the window → no interval.
TEST( HitCollector, WindowBoundaryOnePastL )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1 } );
    tp.add( std::vector<std::uint64_t>{ 7 } ); // 7-1=6 > L=5 → different windows
    HC hc;
    auto r = hc.query( tp, 5, 2 );
    EXPECT_TRUE( r.empty() );
}

// Three lists, M=2: only the stretch where at least two lists overlap is reported.
// list0={1,2,3}, list1={2,3,4}, list2={10,11} (isolated), L=3, M=2
TEST( HitCollector, ThreeListsMPartial )
{
    TP tp( 3 );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3 } );
    tp.add( std::vector<std::uint64_t>{ 2, 3, 4 } );
    tp.add( std::vector<std::uint64_t>{ 10, 11 } ); // isolated; only 1 list → no interval
    HC hc;
    auto r = hc.query( tp, 3, 2 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 4u );
}

// Smoke test: uint32_t instantiation compiles and produces correct results.
TEST( HitCollector, Uint32Instantiation )
{
    TermPostings<std::uint32_t> tp( 2 );
    tp.add( std::vector<std::uint32_t>{ 1, 2, 3 } );
    tp.add( std::vector<std::uint32_t>{ 2, 3, 4 } );
    HitCollector<std::uint32_t, 16> hc;
    auto r = hc.query( tp, 3u, 2 );
    ASSERT_EQ( r.size(), 1u );
    EXPECT_EQ( r[0].left,  1u );
    EXPECT_EQ( r[0].right, 4u );
}

// Out-parameter overload reuses the vector allocation across calls.
TEST( HitCollector, OutParamOverloadReusesVector )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 2, 3 } );
    tp.add( std::vector<std::uint64_t>{ 2, 3, 4 } );
    HC hc;
    std::vector<Interval> out;
    out.reserve( 8 ); // set a capacity that should survive across calls
    hc.query( tp, 3, 2, out );
    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left, 1u );
    EXPECT_GE( out.capacity(), 1u ); // capacity preserved (not shrunk)
}

// A type derived from HitInterval with an extra field, standing in for callers (like
// KmerSeeding::SeedInterval) that need extra per-interval payload beyond what HitCollector
// itself knows about.
struct DerivedInterval : public HC::HitInterval
{
    int payload = -1;
};

// query<IntervalT>() must work generically for a type derived from HitInterval, constructing
// left/right/peak_hits correctly and leaving the derived field at its default.
TEST( HitCollector, GenericQuerySupportsDerivedIntervalType )
{
    TP tp( 2 );
    tp.add( std::vector<std::uint64_t>{ 1, 5 } );
    tp.add( std::vector<std::uint64_t>{ 3, 7 } );
    HC hc;

    std::vector<DerivedInterval> out;
    hc.query( tp, 5, 2, out );
    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,  1u );
    EXPECT_EQ( out[0].right, 7u );
    EXPECT_EQ( out[0].payload, -1 ); // untouched by HitCollector, still at its default

    // The by-value overload cannot deduce IntervalT (it only appears in the return type),
    // so it must be given explicitly here.
    auto const out2 = hc.query<DerivedInterval>( tp, 5, 2 );
    ASSERT_EQ( out2.size(), 1u );
    EXPECT_EQ( out2[0].left,  1u );
    EXPECT_EQ( out2[0].payload, -1 );
}

// =================================================================================================
//     Performance tests
// =================================================================================================

// S=100 lists, 100k entries each, values sparse in [0, 1e9], L=5, M=50.
// Exercises the k-way merge throughput on sparse data.
TEST( HitCollectorPerf, LargeSparseLists )
{
    constexpr int    S        = 100;
    constexpr int    N        = 100'000;
    constexpr std::uint64_t MAX_VAL = 1'000'000'000ULL;
    constexpr std::uint64_t L       = 5;
    constexpr std::size_t   M       = 50;

    std::mt19937_64 rng( 42 );
    std::uniform_int_distribution<std::uint64_t> dist( 0, MAX_VAL );

    TP tp( S );
    for( int s = 0; s < S; ++s ) {
        std::vector<std::uint64_t> v( N );
        for( auto& x : v ) x = dist( rng );
        std::sort( v.begin(), v.end() );
        v.erase( std::unique( v.begin(), v.end() ), v.end() );
        tp.add( std::move( v ) );
    }

    HC hc;
    std::vector<Interval> out;
    auto t0 = std::chrono::high_resolution_clock::now();
    hc.query( tp, L, M, out );
    auto t1 = std::chrono::high_resolution_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>( t1 - t0 ).count();
    std::cout << "[LargeSparseLists] S=" << S << " N~" << N
              << " L=" << L << " M=" << M
              << " → " << out.size() << " intervals, " << ms << " ms\n";
}

// S=10 lists each covering a dense range [0, 99999], L=5, M=5 (all lists).
// Exercises the ring buffer OR-recompute path under maximum hit density.
TEST( HitCollectorPerf, DenseLists )
{
    constexpr int           S = 10;
    constexpr std::uint64_t N = 100'000;
    constexpr std::uint64_t L = 5;
    constexpr std::size_t   M = 5;

    TP tp( S );
    for( int s = 0; s < S; ++s ) {
        std::vector<std::uint64_t> v( N );
        std::iota( v.begin(), v.end(), std::uint64_t{0} );
        tp.add( std::move( v ) );
    }

    HC hc;
    std::vector<Interval> out;
    auto t0 = std::chrono::high_resolution_clock::now();
    hc.query( tp, L, M, out );
    auto t1 = std::chrono::high_resolution_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>( t1 - t0 ).count();
    std::cout << "[DenseLists] S=" << S << " N=" << N
              << " L=" << L << " M=" << M
              << " → " << out.size() << " intervals, " << ms << " ms\n";

    // Dense identical lists: one interval from 0 to N-1.
    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,  0u );
    EXPECT_EQ( out[0].right, N - 1 );
}
