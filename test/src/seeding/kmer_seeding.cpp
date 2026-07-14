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
 * @brief Tests for KmerSeeding.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/inverted_index/builder.hpp"
#include "spear/inverted_index/index.hpp"
#include "spear/seeding/kmer_seeding.hpp"
#include "genesis/util/threading/thread_pool.hpp"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace spear::inverted_index;
using namespace spear::seeding;

// Convenience aliases for the default (uint64_t) instantiation.
using Builder  = InvertedIndexBuilder<>;
using Index    = InvertedIndex<>;
using KS       = KmerSeeding<>;
using Interval = KS::SeedInterval;

// =================================================================================================
//     Shared Fixture Helper
// =================================================================================================

// Builds an index with postings_per_term.size() terms, where term i has the postings given in
// postings_per_term[i] (an empty vector is a valid, present-but-empty term). Writes to `path`.
static void build_index_(
    std::string const& path,
    std::vector<std::vector<std::uint64_t>> const& postings_per_term
) {
    Builder builder( postings_per_term.size() );
    for( std::size_t t = 0; t < postings_per_term.size(); ++t ) {
        for( auto const& pos : postings_per_term[t] ) {
            builder.add( t, pos );
        }
    }
    builder.write( path );
}

// =================================================================================================
//     Construction
// =================================================================================================

TEST( KmerSeeding, ThrowsOnZeroBinWidth )
{
    Index idx; // never opened; the ctor only checks bin_width, not the index
    KS::Config cfg;
    EXPECT_THROW( KS( cfg, idx, 0 ), std::invalid_argument );
}

TEST( KmerSeeding, ConstructsWithNonZeroBinWidth )
{
    Index idx;
    KS::Config cfg;
    EXPECT_NO_THROW( KS( cfg, idx, 100 ) );
}

// =================================================================================================
//     Type Traits
// =================================================================================================

// KmerSeeding holds a reference and owns non-moveable Stats; Query guards per-thread state via
// its destructor. These traits are part of the documented API contract, so pin them down here.
TEST( KmerSeeding, TypeTraits )
{
    EXPECT_FALSE( std::is_copy_constructible_v<KS> );
    EXPECT_FALSE( std::is_move_constructible_v<KS> );
    EXPECT_FALSE( std::is_copy_assignable_v<KS> );
    EXPECT_FALSE( std::is_move_assignable_v<KS> );

    EXPECT_FALSE( std::is_copy_constructible_v<KS::Query> );
    EXPECT_TRUE ( std::is_move_constructible_v<KS::Query> );
    EXPECT_FALSE( std::is_copy_assignable_v<KS::Query> );
    EXPECT_FALSE( std::is_move_assignable_v<KS::Query> );
}

// =================================================================================================
//     Query Lifecycle
// =================================================================================================

TEST( KmerSeeding, StartQueryReentrancyGuard )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }, { 20 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    auto q1 = ks.start_query();
    EXPECT_THROW( (void) ks.start_query(), std::runtime_error );

    std::vector<Interval> out;
    q1.finish_query( 50, out );

    // The guard is released by finish_query(), so starting again must now succeed.
    auto q2 = ks.start_query();
    q2.finish_query( 50, out );
}

TEST( KmerSeeding, QueryMoveConstructTransfersOwnership )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }, { 20 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    auto q1 = ks.start_query();
    q1.add_kmer( 0 );
    auto q2 = std::move( q1 );

    // Moved-from q1 must be a safe no-op on destruction; only q2 needs finishing. If the
    // moved-from destructor were not a no-op, this test would abort the whole process.
    std::vector<Interval> out;
    q2.finish_query( 50, out );
}

// =================================================================================================
//     add_kmer / Truncation
// =================================================================================================

TEST( KmerSeeding, AddKmerAcceptsUpToMaxAndTruncatesBeyond )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10, 20 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    auto q = ks.start_query();
    for( std::size_t i = 0; i < KS::kMaxKmersPerRead; ++i ) {
        EXPECT_TRUE( q.add_kmer( 0 ) ) << "i=" << i;
    }
    EXPECT_FALSE( q.add_kmer( 0 ) );
    EXPECT_FALSE( q.add_kmer( 0 ) ); // still false on repeated over-limit calls

    std::vector<Interval> out;
    q.finish_query( 500, out );

    EXPECT_EQ( ks.stats().truncated_queries, 1u );
}

TEST( KmerSeeding, AddKmerNoTruncationUnderLimit )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }, { 20 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    auto q = ks.start_query();
    EXPECT_TRUE( q.add_kmer( 0 ) );
    EXPECT_TRUE( q.add_kmer( 1 ) );

    std::vector<Interval> out;
    q.finish_query( 50, out );

    EXPECT_EQ( ks.stats().truncated_queries, 0u );
}

// =================================================================================================
//     Window Length / Min Hits Derivation
// =================================================================================================

// bin_width=10, read_length=50 -> window_length = ceil(50/10)+1 = 6, used as the shared fixture
// for the boundary and min-hits tests below.

TEST( KmerSeeding, WindowLengthBoundaryExactlyAtDerivedValue )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 106 }} ); // gap = 6 == window_length
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    std::vector<Interval> out;
    q.finish_query( 50, out );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,  100u );
    EXPECT_EQ( out[0].right, 106u );
}

TEST( KmerSeeding, WindowLengthOnePastDerivedValueYieldsNoInterval )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 107 }} ); // gap = 7 > window_length = 6
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    std::vector<Interval> out;
    q.finish_query( 50, out );

    EXPECT_TRUE( out.empty() );
}

TEST( KmerSeeding, WindowLengthClampedToMax )
{
    // bin_width=1, read_length=1000 -> raw window = ceil(1000/1)+1 = 1001, which must be
    // clamped to kMaxWindowLength. A gap one past kMaxWindowLength must NOT merge under the
    // clamp, even though it trivially would under the unclamped raw value of 1001.
    auto const gap = static_cast<std::uint64_t>( KS::kMaxWindowLength ) + 1;

    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 100 + gap }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 1 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    std::vector<Interval> out;
    q.finish_query( 1000, out );

    EXPECT_TRUE( out.empty() );
}

TEST( KmerSeeding, WindowLengthNotClampedToMax )
{
    // Same as above, but just below the clamp, so the interval should merge.
    auto const gap = static_cast<std::uint64_t>( KS::kMaxWindowLength );

    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 100 + gap }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 1 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    std::vector<Interval> out;
    q.finish_query( 1000, out );

    EXPECT_FALSE( out.empty() );
}

TEST( KmerSeeding, MinHitCountOverridesFraction )
{
    // 3 kmers added; fraction=1.0 would require all 3 to co-occur, but min_hit_count=2 takes
    // priority and only requires 2. Term 2 sits far from 0/1 so it never co-occurs with them;
    // terms 0/1 do co-occur (gap 1 <= window_length 6).
    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 101 }, { 5000 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count    = 2;
    cfg.min_hit_fraction = 1.0; // would require 3 if this were used instead
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    q.add_kmer( 2 );
    std::vector<Interval> out;
    q.finish_query( 50, out );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,       100u );
    EXPECT_EQ( out[0].right,      101u );
    EXPECT_EQ( out[0].peak_hits,  2u );
}

TEST( KmerSeeding, MinHitFractionRounding )
{
    // 5 kmers added, fraction=0.3 -> llround(5*0.3) = llround(1.5) = 2 (round half away from
    // zero), not a truncated 1. With min_hits=2, only the pair at {100,101} forms an interval;
    // the three far-apart singletons never co-occur with anything. If the derivation instead
    // truncated to 1, min_hits=1 would make any non-empty window valid, and the whole span
    // (100..7000) would merge into a single interval with right=7000 instead of 101.
    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 101 }, { 5000 }, { 6000 }, { 7000 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_fraction = 0.3;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    q.add_kmer( 2 );
    q.add_kmer( 3 );
    q.add_kmer( 4 );
    std::vector<Interval> out;
    q.finish_query( 50, out );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,  100u );
    EXPECT_EQ( out[0].right, 101u );
}

TEST( KmerSeeding, MinHitsDerivedFromAddedCountNotFoundCount )
{
    // 4 kmers added: terms 0/1 are found (placed far apart so they never co-occur), terms 2/3
    // are present in the index but have no postings (kEmpty), so they don't consume a
    // TermPostings slot. fraction=0.5 must be computed from count_ (4 added kmers) ->
    // llround(4*0.5)=2, not from list_count() (2 found) -> llround(2*0.5)=1. With min_hits=2
    // and 0/1 never co-occurring, no interval can appear. If the derivation instead used
    // list_count()==2 (-> min_hits=1), the whole span would incorrectly merge into one
    // interval, since min_hits=1 makes any non-empty window valid.
    TempFile tmp;
    build_index_( tmp.path, {{ 100 }, { 5000 }, {}, {}} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_fraction = 0.5;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    q.add_kmer( 2 ); // empty term
    q.add_kmer( 3 ); // empty term
    std::vector<Interval> out;
    q.finish_query( 50, out );

    EXPECT_TRUE( out.empty() );
}

// =================================================================================================
//     Stats
// =================================================================================================

TEST( KmerSeeding, StatsTotalQueriesAccumulates )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    std::vector<Interval> out;
    for( int i = 0; i < 3; ++i ) {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.finish_query( 50, out );
    }

    EXPECT_EQ( ks.stats().total_queries, 3u );
}

TEST( KmerSeeding, StatsEmptyQueriesCountsOnlyEmptyResults )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }, {}} ); // term 0 found, term 1 present-but-empty
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );
    std::vector<Interval> out;

    // Non-empty result: single found list, default min_hits=1.
    {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.finish_query( 50, out );
    }
    EXPECT_FALSE( out.empty() );

    // Empty result: only the empty term is queried, so list_count()==0.
    {
        auto q = ks.start_query();
        q.add_kmer( 1 );
        q.finish_query( 50, out );
    }
    EXPECT_TRUE( out.empty() );

    EXPECT_EQ( ks.stats().total_queries, 2u );
    EXPECT_EQ( ks.stats().empty_queries, 1u );
}

TEST( KmerSeeding, StatsTruncatedQueriesCountsOncePerQuery )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );
    std::vector<Interval> out;

    // Truncated query: exceed kMaxKmersPerRead twice; must still only count once.
    {
        auto q = ks.start_query();
        for( std::size_t i = 0; i < KS::kMaxKmersPerRead; ++i ) {
            q.add_kmer( 0 );
        }
        q.add_kmer( 0 );
        q.add_kmer( 0 );
        q.finish_query( 500, out );
    }

    // Non-truncated query.
    {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.finish_query( 50, out );
    }

    EXPECT_EQ( ks.stats().truncated_queries, 1u );
}

TEST( KmerSeeding, StatsSoftCappedKmersAccumulatesAcrossQueries )
{
    // Term 0 has 3 postings; a soft cap of 1 makes every add_kmer(0) count as soft-capped.
    TempFile tmp;
    build_index_( tmp.path, {{ 10, 20, 30 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.max_occurrences_per_kmer = 1;
    KS ks( cfg, idx, 100 );

    std::vector<Interval> out;
    {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.finish_query( 50, out );
    }
    {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.add_kmer( 0 );
        q.finish_query( 50, out );
    }

    EXPECT_EQ( ks.stats().soft_capped_kmers, 3u );
}

TEST( KmerSeeding, StatsPostingsRollupCountsFoundAbsentAndHardCapped )
{
    // term 0: 3 postings, hard-capped at build time (max_postings_per_term=1).
    // term 1: 1 posting, found normally.
    // term 2: present but empty (kEmpty).
    TempFile tmp;
    {
        Builder::Config cfg;
        cfg.max_postings_per_term = 1;
        Builder builder( 3, cfg );
        builder.add( 0, 10 );
        builder.add( 0, 20 );
        builder.add( 0, 30 );
        builder.add( 1, 40 );
        builder.write( tmp.path );
    }
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );
    std::vector<Interval> out;

    auto q = ks.start_query();
    q.add_kmer( 0 ); // hard-capped
    q.add_kmer( 1 ); // found
    q.add_kmer( 2 ); // absent (empty)
    q.finish_query( 50, out );

    auto const s = ks.stats();
    EXPECT_EQ( s.found_kmers,       1u );
    EXPECT_EQ( s.absent_kmers,      1u );
    EXPECT_EQ( s.hard_capped_kmers, 1u );
}

TEST( KmerSeeding, StatsNumKmersHistogramTracksAddedCount )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10 }, { 20 }, { 30 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );
    std::vector<Interval> out;

    { // 0 k-mers added
        auto q = ks.start_query();
        q.finish_query( 50, out );
    }
    { // 2 k-mers added
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.add_kmer( 1 );
        q.finish_query( 50, out );
    }
    { // 3 k-mers added
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.add_kmer( 1 );
        q.add_kmer( 2 );
        q.finish_query( 50, out );
    }

    auto const hist = ks.stats().num_kmers_histogram;
    ASSERT_EQ( hist.size(), KS::kMaxKmersPerRead + 1 );
    EXPECT_EQ( hist[0], 1u );
    EXPECT_EQ( hist[2], 1u );
    EXPECT_EQ( hist[3], 1u );
}

TEST( KmerSeeding, StatsNumSeedsHistogramTracksIntervalCount )
{
    // Three far-apart pairs, each forming exactly one interval under min_hit_count=2 and
    // window_length=5 (bin_width=10, read_length=40 -> ceil(40/10)+1=5) -- mirrors the
    // concurrency fixture's profile construction (build_concurrent_profiles_ below). All three
    // pairs are queried together, so out should contain exactly 3 intervals: large enough to
    // prove the histogram actually grows past its initial (empty) size, not just that bin 0/1
    // work as edge cases.
    TempFile tmp;
    Builder builder( 6 );
    for( std::size_t p = 0; p < 3; ++p ) {
        auto const base = p * 1000;
        builder.add( static_cast<Index::term_index_type>( 2 * p ),     base + 1 );
        builder.add( static_cast<Index::term_index_type>( 2 * p ),     base + 5 );
        builder.add( static_cast<Index::term_index_type>( 2 * p + 1 ), base + 3 );
        builder.add( static_cast<Index::term_index_type>( 2 * p + 1 ), base + 7 );
    }
    builder.write( tmp.path );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    for( Index::term_index_type t = 0; t < 6; ++t ) {
        q.add_kmer( t );
    }
    std::vector<Interval> out;
    q.finish_query( 40, out );
    ASSERT_EQ( out.size(), 3u );

    auto const hist = ks.stats().num_seeds_histogram;
    ASSERT_GT( hist.size(), 3u );
    EXPECT_EQ( hist[3], 1u );
}

TEST( KmerSeeding, StatsTopPeakHitsHistogramBin0ForEmptyResultAndPeakForNonEmpty )
{
    // term 2 alone gives list_count()==1 < min_hit_count==2 -> empty result -> bin 0.
    // terms 0/1 reuse the EndToEnd fixture (list0={1,5}, list1={3,7}) -> peak_hits==2 -> bin 2.
    TempFile tmp;
    build_index_( tmp.path, {{ 1, 5 }, { 3, 7 }, { 999999 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );
    std::vector<Interval> out;

    {
        auto q = ks.start_query();
        q.add_kmer( 2 );
        q.finish_query( 40, out );
    }
    EXPECT_TRUE( out.empty() );

    {
        auto q = ks.start_query();
        q.add_kmer( 0 );
        q.add_kmer( 1 );
        q.finish_query( 40, out );
    }
    ASSERT_FALSE( out.empty() );
    ASSERT_EQ( out[0].peak_hits, 2u );

    auto const hist = ks.stats().top_peak_hits_histogram;
    EXPECT_EQ( hist[0], 1u );
    EXPECT_EQ( hist[2], 1u );
}

TEST( KmerSeeding, CollectStatsFalseSkipsBookkeepingButNotSeeding )
{
    // Reuse the EndToEnd fixture (list0={1,5}, list1={3,7}) so the seeding result itself
    // (out) has a known expected value to check against, proving collect_stats=false only
    // switches off the stats bookkeeping and doesn't affect the actual seeding logic.
    TempFile tmp;
    build_index_( tmp.path, {{ 1, 5 }, { 3, 7 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count   = 2;
    cfg.collect_stats   = false;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    std::vector<Interval> out;
    q.finish_query( 40, out );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,      1u );
    EXPECT_EQ( out[0].right,     7u );
    EXPECT_EQ( out[0].peak_hits, 2u );

    auto const s = ks.stats();
    EXPECT_EQ( s.total_queries, 0u );
    EXPECT_EQ( s.found_kmers,   0u );
    EXPECT_TRUE( s.num_seeds_histogram.empty() );
}

// =================================================================================================
//     End-to-End Functional
// =================================================================================================

TEST( KmerSeeding, EndToEndMatchesHitCollectorReference )
{
    // Mirrors HitCollector's TwoListsOneInterval fixture (list0={1,5}, list1={3,7}, L=5, M=2)
    // through the full KmerSeeding API: bin_width=10, read_length=40 -> window_length =
    // ceil(40/10)+1 = 5.
    TempFile tmp;
    build_index_( tmp.path, {{ 1, 5 }, { 3, 7 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    auto out = q.finish_query( 40 ); // by-value overload

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,      1u );
    EXPECT_EQ( out[0].right,     7u );
    EXPECT_EQ( out[0].peak_hits, 2u );
}

TEST( KmerSeeding, FinishQueryOutParamClearsPriorContents )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 1, 5 }, { 3, 7 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    std::vector<Interval> out;
    out.push_back( Interval{ 999, 999, 999 } ); // stale contents from a prior round

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    q.finish_query( 40, out );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left, 1u );
}

// =================================================================================================
//     RunQuery Convenience Overload
// =================================================================================================

TEST( KmerSeeding, RunQueryMatchesManualQuery )
{
    // Same fixture as EndToEndMatchesHitCollectorReference: run_query() must produce identical
    // output to driving start_query()/add_kmer()/finish_query() manually.
    TempFile tmp;
    build_index_( tmp.path, {{ 1, 5 }, { 3, 7 }} );
    Index idx;
    idx.open( tmp.path );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, 10 );

    std::vector<KS::kmer_index_type> const kmers{ 0, 1 };
    auto const out = ks.run_query( kmers, 40 ); // by-value overload

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,      1u );
    EXPECT_EQ( out[0].right,     7u );
    EXPECT_EQ( out[0].peak_hits, 2u );
}

TEST( KmerSeeding, RunQueryTruncatesBeyondMaxKmersPerRead )
{
    TempFile tmp;
    build_index_( tmp.path, {{ 10, 20 }} );
    Index idx;
    idx.open( tmp.path );

    KS ks( {}, idx, 100 );

    std::vector<KS::kmer_index_type> kmers( KS::kMaxKmersPerRead + 2, 0 );
    std::vector<Interval> out;
    ks.run_query( kmers, 500, out );

    EXPECT_EQ( ks.stats().truncated_queries, 1u );
}

// =================================================================================================
//     PositionT Instantiation
// =================================================================================================

TEST( KmerSeeding, Uint32Instantiation )
{
    using Builder32 = InvertedIndexBuilder<std::uint32_t>;
    using Index32   = InvertedIndex<std::uint32_t>;
    using KS32      = KmerSeeding<std::uint32_t>;

    TempFile tmp;
    {
        Builder32 builder( 2 );
        builder.add( 0, 1 );
        builder.add( 0, 5 );
        builder.add( 1, 3 );
        builder.add( 1, 7 );
        builder.write( tmp.path );
    }
    Index32 idx;
    idx.open( tmp.path );

    KS32::Config cfg;
    cfg.min_hit_count = 2;
    KS32 ks( cfg, idx, 10 );

    auto q = ks.start_query();
    q.add_kmer( 0 );
    q.add_kmer( 1 );
    auto out = q.finish_query( 40 );

    ASSERT_EQ( out.size(), 1u );
    EXPECT_EQ( out[0].left,  1u );
    EXPECT_EQ( out[0].right, 7u );
}

// =================================================================================================
//     Concurrency
// =================================================================================================

// A profile is a pair of term indices whose postings are crafted to produce a single, exactly
// known SeedInterval under kConcurrentBinWidth/kConcurrentReadLen (window_length = 5) and
// min_hit_count = 2 -- mirrors HitCollector's TwoListsOneInterval fixture, offset by
// kProfileSpacing per profile so that different profiles' postings are far apart and never
// interact.
namespace {

struct Profile
{
    Index::term_index_type term_a;
    Index::term_index_type term_b;
    std::uint64_t           expected_left;
    std::uint64_t           expected_right;
};

} // anonymous namespace

static constexpr std::uint64_t kProfileSpacing     = 1000;
static constexpr std::size_t   kConcurrentBinWidth = 10;
static constexpr std::size_t   kConcurrentReadLen  = 40; // -> window_length = ceil(40/10)+1 = 5

static std::vector<Profile> build_concurrent_profiles_( Builder& builder, std::size_t num_profiles )
{
    std::vector<Profile> profiles;
    profiles.reserve( num_profiles );
    for( std::size_t p = 0; p < num_profiles; ++p ) {
        auto const base   = static_cast<std::uint64_t>( p ) * kProfileSpacing;
        auto const term_a = 2 * p;
        auto const term_b = 2 * p + 1;
        builder.add( term_a, base + 1 );
        builder.add( term_a, base + 5 );
        builder.add( term_b, base + 3 );
        builder.add( term_b, base + 7 );
        profiles.push_back({ term_a, term_b, base + 1, base + 7 });
    }
    return profiles;
}

// Every thread repeatedly runs the same fixed profile against a shared engine and asserts the
// exact expected result on every round. A leak of thread_local state between threads (the
// per-thread TermPostings/HitCollector in start_query()) would show up as a sporadically wrong
// result for some profile, not merely as a crash -- which is what makes this stronger than a
// "no crash" style stress test.
static void run_concurrent_deterministic_( std::size_t num_threads, std::size_t iterations )
{
    using namespace genesis::util::threading;

    constexpr std::size_t kNumProfiles = 4;

    TempFile tmp;
    std::vector<Profile> profiles;
    {
        Builder builder( 2 * kNumProfiles );
        profiles = build_concurrent_profiles_( builder, kNumProfiles );
        builder.write( tmp.path );
    }

    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, kConcurrentBinWidth );

    auto pool = std::make_shared<ThreadPool>( num_threads );
    std::vector<ProactiveFuture<bool>> futures;
    futures.reserve( num_threads );

    for( std::size_t t = 0; t < num_threads; ++t ) {
        auto const& profile = profiles[ t % kNumProfiles ];
        futures.emplace_back( pool->enqueue_and_retrieve(
            [&ks, &profile, iterations]() -> bool
            {
                std::vector<Interval> out;
                for( std::size_t i = 0; i < iterations; ++i ) {
                    auto q = ks.start_query();
                    q.add_kmer( profile.term_a );
                    q.add_kmer( profile.term_b );
                    q.finish_query( kConcurrentReadLen, out );
                    if(
                        out.size() != 1 ||
                        out[0].left != profile.expected_left ||
                        out[0].right != profile.expected_right ||
                        out[0].peak_hits != 2
                    ) {
                        return false;
                    }
                }
                return true;
            }
        ));
    }

    for( auto& f : futures ) {
        EXPECT_TRUE( f.get() ) << "num_threads=" << num_threads << " iterations=" << iterations;
    }

    EXPECT_EQ( ks.stats().total_queries, num_threads * iterations );
    EXPECT_EQ( ks.stats().empty_queries, 0u );
    EXPECT_EQ( ks.stats().truncated_queries, 0u );

    // Every query yields exactly one interval, so all concurrent increments land on the same
    // bin (index 1) of the mutex-guarded num_seeds_histogram_ -- this is what actually exercises
    // that mechanism under real concurrent growth/increment, not just the single-threaded logic
    // covered by the StatsNumSeedsHistogram* tests above.
    auto const seeds_hist = ks.stats().num_seeds_histogram;
    ASSERT_EQ( seeds_hist.size(), 2u );
    EXPECT_EQ( seeds_hist[0], 0u );
    EXPECT_EQ( seeds_hist[1], num_threads * iterations );
}

// Sweeps thread counts (mirroring InvertedIndex's StressNoCap/StressCapped convention) and
// repeats several times, since races that corrupt results are often scheduling-order dependent
// and may not manifest on every run.
TEST( KmerSeedingConcurrent, Deterministic )
{
    constexpr std::size_t kIterations = 75;
    constexpr std::size_t kRepeats    = 4;
    for( std::size_t r = 0; r < kRepeats; ++r ) {
        for( std::size_t num_threads : { 1u, 2u, 4u, 8u } ) {
            run_concurrent_deterministic_( num_threads, kIterations );
        }
    }
}

// Large-scale concurrent soak test, mirroring InvertedIndex's DISABLED_StressPerf: a big
// synthetic index, sustained concurrent query load, and throughput reported to stdout. Disabled
// by default; enable with --gtest_also_run_disabled_tests. Raise kNumProfiles / kNumQueries for
// a longer soak; lower for a quick smoke run. Correctness here is light-touch (total query count
// only) since ConcurrentDeterministic above already covers exact per-query correctness; this
// test exists to observe throughput and contention under sustained load.
//
// One task is enqueued per query, rather than a fixed small number of long-running worker
// loops: with only as many tasks as pool threads, each worker just runs one big closure to
// completion and the pool's queue is barely touched, which says little about the threading
// model itself. Submitting one detached task per query instead exercises the pool's actual
// production usage pattern -- a single producer thread enqueuing one task per unit of work
// (mirroring commands/map/align.cpp, which enqueues one detached task per read), consumed by
// kNumThreads workers pulling from the shared queue -- and puts real, sustained pressure on the
// queue's enqueue/dequeue synchronization instead of just running kNumThreads independent loops.
// TEST( KmerSeedingConcurrent, StressPerf )
TEST( KmerSeedingConcurrent, DISABLED_StressPerf )
{
    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::duration<double, std::milli>;
    using namespace genesis::util::threading;

    constexpr std::size_t kNumProfiles = 2'500'000;   // -> 2 * kNumProfiles terms
    constexpr std::size_t kNumThreads  = 8;           // pool worker count
    constexpr std::size_t kNumQueries  = 100'000'000; // one task enqueued per query

    std::cout << "\n[KmerSeedingConcurrentStressPerf] profiles=" << kNumProfiles
              << "  terms=" << 2 * kNumProfiles
              << "  threads=" << kNumThreads
              << "  queries=" << kNumQueries << "\n";

    TempFile tmp;
    std::vector<Profile> profiles;
    auto const t_build0 = Clock::now();
    {
        Builder builder( 2 * kNumProfiles );
        profiles = build_concurrent_profiles_( builder, kNumProfiles );
        builder.write( tmp.path );
    }
    double const build_ms = Ms( Clock::now() - t_build0 ).count();
    std::cout << "[KmerSeedingConcurrentStressPerf] build+write: " << build_ms / 1e3 << " s\n";

    auto const t_precompute0 = Clock::now();
    Index idx;
    idx.open( tmp.path, Index::OpenMode::kLoadAll );

    KS::Config cfg;
    cfg.min_hit_count = 2;
    KS ks( cfg, idx, kConcurrentBinWidth );

    // Precompute which profile each query uses with a single RNG pass up front, so that each
    // enqueued task is a small, self-contained unit of work rather than needing its own RNG.
    // Also precompute the checksum we expect the queries to produce (see result_checksum below),
    // since every profile query is designed to yield exactly one interval with a known
    // left/right/peak_hits.
    std::vector<std::size_t> profile_indices( kNumQueries );
    std::uint64_t expected_checksum = 0;
    {
        std::mt19937_64 rng( 42 );
        std::uniform_int_distribution<std::size_t> dist( 0, profiles.size() - 1 );
        for( auto& pi : profile_indices ) {
            pi = dist( rng );
            auto const& profile = profiles[pi];
            expected_checksum += profile.expected_left + profile.expected_right + 2; // peak_hits
        }
    }
    double const precompute_ms = Ms( Clock::now() - t_precompute0 ).count();
    std::cout << "[KmerSeedingConcurrentStressPerf] precompute: " << precompute_ms / 1e3 << " s\n";

    auto pool = std::make_shared<ThreadPool>( kNumThreads );

    // enqueue_detached() cannot propagate exceptions back to the caller; without this callback,
    // a bug that makes a query throw would terminate the whole test process instead of failing
    // the assertion below.
    std::atomic<bool> task_failed{ false };
    pool->detached_task_exception_callback = [&task_failed]( std::exception_ptr ) {
        task_failed = true;
    };

    // Accumulated from every task's actual result below. Without this, `out` is written by
    // finish_query() but never read by anything, so an optimizer -- especially under this
    // build's LTO/IPO -- could in principle prove the query work has no observable effect and
    // eliminate it. Folding the results into an atomic that feeds the checksum printed (and
    // asserted against the precomputed expected_checksum) below makes the work observable, and
    // as a bonus turns this from a "did it crash" soak test into one that also verifies the
    // seeding logic still returns correct results under sustained concurrent load.
    std::atomic<std::uint64_t> result_checksum{ 0 };

    auto const t_query0 = Clock::now();
    for( std::size_t i = 0; i < kNumQueries; ++i ) {
        pool->enqueue_detached(
            [&ks, &profiles, &profile_indices, &result_checksum, i]()
            {
                auto const& profile = profiles[ profile_indices[i] ];
                thread_local std::vector<Interval> out;
                auto q = ks.start_query();
                q.add_kmer( profile.term_a );
                q.add_kmer( profile.term_b );
                q.finish_query( kConcurrentReadLen, out );

                std::uint64_t local_sum = 0;
                for( auto const& iv : out ) {
                    local_sum += iv.left + iv.right + iv.peak_hits;
                }
                result_checksum.fetch_add( local_sum, std::memory_order_relaxed );
            }
        );
    }
    pool->wait_for_all_pending_tasks();
    double const query_ms = Ms( Clock::now() - t_query0 ).count();

    std::cout << "[KmerSeedingConcurrentStressPerf] query: " << query_ms / 1e3 << " s"
              << "  (" << kNumQueries / ( query_ms / 1e3 ) / 1e6 << " M queries/s)"
              << "  checksum=" << result_checksum.load() << "\n";

    EXPECT_FALSE( task_failed.load() );
    EXPECT_EQ( ks.stats().total_queries, kNumQueries );
    EXPECT_EQ( result_checksum.load(), expected_checksum );
}
