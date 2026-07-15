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
 * @brief Tests for filter_seed_intervals() / SeedFilterConfig.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/seeding/seed_filter.hpp"

#include <cstddef>
#include <vector>

using namespace spear::seeding;

// =================================================================================================
//     Test fixture type
// =================================================================================================

// Minimal stand-in for KmerSeeding<PositionT>::SeedInterval: filter_seed_intervals() is a
// template that only ever touches peak_hits, so a tiny local type both keeps this test free of
// any InvertedIndex/HitCollector setup and demonstrates the function's genericity.
struct TestInterval
{
    int         left       = 0;
    int         right      = 0;
    std::size_t peak_hits  = 0;
};

namespace {

std::vector<TestInterval> make_intervals( std::vector<std::size_t> const& peak_hits )
{
    std::vector<TestInterval> result;
    result.reserve( peak_hits.size() );
    for( auto const ph : peak_hits ) {
        result.push_back( TestInterval{ 0, 0, ph } );
    }
    return result;
}

std::vector<std::size_t> peak_hits_of( std::vector<TestInterval> const& intervals )
{
    std::vector<std::size_t> result;
    result.reserve( intervals.size() );
    for( auto const& iv : intervals ) {
        result.push_back( iv.peak_hits );
    }
    return result;
}

} // anonymous namespace

// =================================================================================================
//     Basic behavior
// =================================================================================================

TEST( SeedFilter, EmptyInputStaysEmpty )
{
    std::vector<TestInterval> intervals;
    filter_seed_intervals( intervals, SeedFilterConfig{} );
    EXPECT_TRUE( intervals.empty() );
}

TEST( SeedFilter, NoFiltersActiveKeepsEverything )
{
    auto intervals = make_intervals( { 10, 8, 5, 1 } );
    filter_seed_intervals( intervals, SeedFilterConfig{} );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8, 5, 1 } ));
}

// =================================================================================================
//     max_seeds
// =================================================================================================

TEST( SeedFilter, MaxSeedsCapsCount )
{
    auto intervals = make_intervals( { 10, 8, 5, 1 } );
    SeedFilterConfig cfg;
    cfg.max_seeds = 2;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8 } ));
}

TEST( SeedFilter, MaxSeedsLargerThanSizeIsNoOp )
{
    auto intervals = make_intervals( { 10, 8 } );
    SeedFilterConfig cfg;
    cfg.max_seeds = 100;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8 } ));
}

TEST( SeedFilter, MaxSeedsZeroIsDisabled )
{
    auto intervals = make_intervals( { 10, 8, 5 } );
    SeedFilterConfig cfg;
    cfg.max_seeds = 0;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8, 5 } ));
}

// =================================================================================================
//     min_peak_hits
// =================================================================================================

TEST( SeedFilter, MinPeakHitsDropsBelowFloor )
{
    auto intervals = make_intervals( { 10, 8, 5, 1 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits = 5;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8, 5 } ));
}

TEST( SeedFilter, MinPeakHitsExactBoundaryIsKept )
{
    // peak_hits == min_peak_hits must survive (only strictly-below values are dropped).
    auto intervals = make_intervals( { 10, 5 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits = 5;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 5 } ));
}

TEST( SeedFilter, MinPeakHitsZeroIsDisabled )
{
    auto intervals = make_intervals( { 10, 1 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits = 0;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 1 } ));
}

// =================================================================================================
//     peak_hits_window
// =================================================================================================

TEST( SeedFilter, PeakHitsWindowDropsTooFarBelowBest )
{
    auto intervals = make_intervals( { 10, 8, 5, 1 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_window = 3;  // keep peak_hits >= 10 - 3 = 7
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8 } ));
}

TEST( SeedFilter, PeakHitsWindowExactBoundaryIsKept )
{
    // best - peak_hits == window must survive (only strictly-greater gaps are dropped).
    auto intervals = make_intervals( { 10, 7, 6 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_window = 3;  // keep peak_hits >= 7
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 7 } ));
}

TEST( SeedFilter, PeakHitsWindowZeroKeepsOnlyTiesWithBest )
{
    auto intervals = make_intervals( { 10, 10, 9 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_window = 0;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 10 } ));
}

TEST( SeedFilter, PeakHitsWindowNegativeIsDisabled )
{
    auto intervals = make_intervals( { 10, 1 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_window = -1;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 1 } ));
}

// =================================================================================================
//     peak_hits_fraction
// =================================================================================================

TEST( SeedFilter, PeakHitsFractionDropsBelowRatio )
{
    auto intervals = make_intervals( { 10, 8, 4, 1 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_fraction = 0.5;  // keep peak_hits >= 5
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8 } ));
}

TEST( SeedFilter, PeakHitsFractionExactBoundaryIsKept )
{
    // peak_hits == fraction * best must survive (only strictly-below values are dropped).
    auto intervals = make_intervals( { 10, 5 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_fraction = 0.5;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 5 } ));
}

TEST( SeedFilter, PeakHitsFractionZeroIsDisabled )
{
    auto intervals = make_intervals( { 10, 1 } );
    SeedFilterConfig cfg;
    cfg.peak_hits_fraction = 0.0;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 1 } ));
}

// =================================================================================================
//     Combinations
// =================================================================================================

TEST( SeedFilter, MultipleFiltersCombineAsAnd )
{
    // best=10; peak_hits_window=4 keeps >=6 (10,8,6); max_seeds=2 further caps to the first 2.
    auto intervals = make_intervals( { 10, 8, 6, 3 } );
    SeedFilterConfig cfg;
    cfg.max_seeds         = 2;
    cfg.peak_hits_window  = 4;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8 } ));
}

TEST( SeedFilter, MostRestrictiveFilterWinsRegardlessOfOrder )
{
    // min_peak_hits alone would keep {10,8,6,4}; peak_hits_fraction alone would keep {10,8,6}
    // (>= 0.5*10 = 5); the combination must match the more restrictive one.
    auto intervals = make_intervals( { 10, 8, 6, 4, 2 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits      = 3;
    cfg.peak_hits_fraction = 0.5;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 10, 8, 6 } ));
}

TEST( SeedFilter, SingleElementSurvivesWhenPassing )
{
    auto intervals = make_intervals( { 5 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits = 5;
    filter_seed_intervals( intervals, cfg );
    EXPECT_EQ( peak_hits_of( intervals ), ( std::vector<std::size_t>{ 5 } ));
}

TEST( SeedFilter, EverythingFilteredOutLeavesEmptyVector )
{
    auto intervals = make_intervals( { 10, 8, 5 } );
    SeedFilterConfig cfg;
    cfg.min_peak_hits = 100;
    filter_seed_intervals( intervals, cfg );
    EXPECT_TRUE( intervals.empty() );
}

// =================================================================================================
//     merge_seed_intervals
// =================================================================================================

TEST( SeedMerge, BothEmptyStaysEmpty )
{
    std::vector<TestInterval> a;
    std::vector<TestInterval> b;
    auto const merged = merge_seed_intervals( a, b );
    EXPECT_TRUE( merged.empty() );
}

TEST( SeedMerge, OneEmptyPassesOtherThrough )
{
    auto a = make_intervals( { 10, 5, 1 } );
    std::vector<TestInterval> b;
    auto const merged = merge_seed_intervals( a, b );
    EXPECT_EQ( peak_hits_of( merged ), ( std::vector<std::size_t>{ 10, 5, 1 } ));

    auto const merged2 = merge_seed_intervals( b, a );
    EXPECT_EQ( peak_hits_of( merged2 ), ( std::vector<std::size_t>{ 10, 5, 1 } ));
}

TEST( SeedMerge, InterleavesWhilePreservingDescendingOrder )
{
    auto a = make_intervals( { 9, 6, 3 } );
    auto b = make_intervals( { 10, 5, 1 } );
    auto const merged = merge_seed_intervals( a, b );
    EXPECT_EQ( peak_hits_of( merged ), ( std::vector<std::size_t>{ 10, 9, 6, 5, 3, 1 } ));
}

TEST( SeedMerge, ResultIsValidPrefixInputForFilterSeedIntervals )
{
    // Confirms the merge-then-filter pattern used by align.cpp: the merged list stays sorted
    // descending, which filter_seed_intervals() requires (it relies on intervals.front() being
    // the maximum, and on a prefix scan).
    auto a = make_intervals( { 8, 4 } );
    auto b = make_intervals( { 10, 6, 2 } );
    auto merged = merge_seed_intervals( a, b );

    SeedFilterConfig cfg;
    cfg.max_seeds = 3;
    filter_seed_intervals( merged, cfg );
    EXPECT_EQ( peak_hits_of( merged ), ( std::vector<std::size_t>{ 10, 8, 6 } ));
}

TEST( SeedMerge, OutParamOverloadClearsPriorContents )
{
    auto a = make_intervals( { 5, 2 } );
    auto b = make_intervals( { 4, 1 } );
    std::vector<TestInterval> out = make_intervals( { 999 } ); // stale contents from a prior round
    merge_seed_intervals( a, b, out );
    EXPECT_EQ( peak_hits_of( out ), ( std::vector<std::size_t>{ 5, 4, 2, 1 } ));
}
