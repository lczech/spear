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

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "spear/inverted_index/builder.hpp"

using namespace spear::inverted_index;

// We just test the default int sizes here. The template parameters are mostly
// for use case optimization and do not affect the public API and behavior.
using Builder = InvertedIndexBuilder<>;

// =================================================================================================
//     Basic Tests
// =================================================================================================

TEST(InvertedIndexBuilder, ConstructorValidationAndBoundsChecks)
{
    {
        Builder::Config cfg;
        cfg.pending_capacity = static_cast<std::size_t>(
            std::numeric_limits<Builder::pending_count_type>::max()
        ) + 1;
        EXPECT_THROW(Builder(4, cfg), std::invalid_argument);
    }

    {
        Builder::Config cfg;
        cfg.max_postings_per_term = std::numeric_limits<Builder::stored_count_type>::max();
        EXPECT_THROW(Builder(4, cfg), std::invalid_argument);
    }

    {
        Builder builder(3);
        EXPECT_EQ(builder.num_term_indices(), 3);
        EXPECT_THROW((void) builder.add(3, 42), std::out_of_range);
        EXPECT_THROW((void) builder.postings(3), std::out_of_range);
        EXPECT_THROW((void) builder.is_capped(3), std::out_of_range);
    }
}

TEST(InvertedIndexBuilder, FinalizeMergesSortsDeduplicatesAcrossFlushesAndTerms)
{
    Builder::Config cfg;
    cfg.pending_capacity = 3;
    cfg.max_postings_per_term = 100;

    Builder builder(4, cfg);

    // Term 0:
    // - adjacent duplicate in pending buffer ignored immediately
    // - overflow flush includes triggering element
    // - later finalize merges pending with already-compressed postings
    builder.add(0, 50);
    builder.add(0, 50);
    builder.add(0, 10);
    builder.add(0, 30);
    builder.add(0, 20);  // overflow-triggered flush of {50,10,30,20}
    builder.add(0, 30);
    builder.add(0, 40);

    // Term 1:
    // - separate scope / independent entry
    // - unsorted inserts and duplicates across flush boundary
    builder.add(1, 7);
    builder.add(1, 3);
    builder.add(1, 7);
    builder.add(1, 1);   // flush
    builder.add(1, 5);
    builder.add(1, 3);

    builder.finalize();

    EXPECT_FALSE(builder.is_capped(0));
    EXPECT_FALSE(builder.is_capped(1));
    EXPECT_FALSE(builder.is_capped(2));

    EXPECT_EQ(builder.posting_count(0), 5);
    EXPECT_EQ(builder.posting_count(1), 4);
    EXPECT_EQ(builder.posting_count(2), 0);

    EXPECT_EQ(builder.postings(0), (std::vector<Builder::position_type>{10, 20, 30, 40, 50}));
    EXPECT_EQ(builder.postings(1), (std::vector<Builder::position_type>{1, 3, 5, 7}));
    EXPECT_EQ(builder.postings(2), (std::vector<Builder::position_type>{}));
}

TEST(InvertedIndexBuilder, CappingDiscardsPostingsAndIgnoresFutureAdds)
{
    Builder::Config cfg;
    cfg.pending_capacity = 2;
    cfg.max_postings_per_term = 4;

    Builder builder(2, cfg);

    builder.add(0, 10);
    builder.add(0, 20);
    builder.add(0, 30);  // flush -> {10,20,30}

    builder.add(0, 40);
    builder.add(0, 50);  // finalize will exceed cap -> capped

    builder.finalize();

    EXPECT_TRUE(builder.is_capped(0));
    EXPECT_EQ(builder.posting_count(0), builder.capped_sentinel());
    EXPECT_EQ(builder.postings(0), (std::vector<Builder::position_type>{}));

    // Once capped, future adds must be ignored permanently.
    builder.add(0, 60);
    builder.add(0, 70);
    builder.finalize();

    EXPECT_TRUE(builder.is_capped(0));
    EXPECT_EQ(builder.posting_count(0), builder.capped_sentinel());
    EXPECT_EQ(builder.postings(0), (std::vector<Builder::position_type>{}));

    // Other entries still behave normally.
    builder.add(1, 9);
    builder.add(1, 2);
    builder.add(1, 9);
    builder.finalize();

    EXPECT_FALSE(builder.is_capped(1));
    EXPECT_EQ(builder.posting_count(1), 2);
    EXPECT_EQ(builder.postings(1), (std::vector<Builder::position_type>{2, 9}));
}

// =================================================================================================
//     Performance Tests
// =================================================================================================

TEST(InvertedIndexBuilder, PendingCapacityPerformance)
{
    using clock = std::chrono::steady_clock;
    using ms_t  = std::chrono::duration<double, std::milli>;

    // Keep this reasonably large so timing differences show up,
    // but not so large that the test becomes annoying to run.
    constexpr std::size_t kNumTerms = 128;
    constexpr std::size_t kPerTerm  = 512;

    std::vector<std::size_t> const pending_capacities = {0, 1, 2, 4, 8, 16};

    for( auto const pending_capacity : pending_capacities ) {
        Builder::Config cfg;
        cfg.pending_capacity       = pending_capacity;
        cfg.max_postings_per_term  = kPerTerm + 16; // comfortably above inserted count

        Builder builder(kNumTerms, cfg);

        auto const start = clock::now();

        // Deterministic, strictly monotonic positions per term.
        // For each term t, the inserted positions are:
        //   t * stride + 0, t * stride + 1, ..., t * stride + (kPerTerm - 1)
        // This is strictly increasing within each term, unique across terms,
        // and cheap to compute.
        constexpr std::uint64_t stride = 1ULL << 32;

        for( std::size_t i = 0; i < kPerTerm; ++i ) {
            for( std::size_t t = 0; t < kNumTerms; ++t ) {
                auto const pos = static_cast<Builder::position_type>(t * stride + i);
                builder.add(static_cast<Builder::term_index_type>(t), pos);
            }
        }

        builder.finalize();

        auto const end = clock::now();
        auto const elapsed_ms = std::chrono::duration_cast<ms_t>(end - start).count();

        LOG_MSG << "pending_capacity=" << pending_capacity
                << " time=" << elapsed_ms << " ms";

        // Correctness check afterwards: every term must contain exactly the
        // inserted monotonic sequence.
        for( std::size_t t = 0; t < kNumTerms; ++t ) {
            EXPECT_FALSE(builder.is_capped(static_cast<Builder::term_index_type>(t)));
            EXPECT_EQ(
                builder.posting_count(static_cast<Builder::term_index_type>(t)),
                kPerTerm
            );

            auto const postings =
                builder.postings(static_cast<Builder::term_index_type>(t));

            ASSERT_EQ(postings.size(), kPerTerm);

            for( std::size_t i = 0; i < kPerTerm; ++i ) {
                auto const expected =
                    static_cast<Builder::position_type>(t * stride + i);
                EXPECT_EQ(postings[i], expected)
                    << "pending_capacity=" << pending_capacity
                    << ", term=" << t
                    << ", index=" << i;
            }
        }
    }
}
