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
 * @brief Tests for ReferenceCollection and its free functions.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/reference/collection.hpp"
#include "spear/reference/functions.hpp"
#include "genesis/sequence/sequence.hpp"

#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace spear::reference;

// =================================================================================================
//     Construction and Accessors
// =================================================================================================

TEST(ReferenceCollection, EmptyCollection)
{
    ReferenceCollection col;
    EXPECT_TRUE( col.empty() );
    EXPECT_EQ( col.size(), 0u );
    EXPECT_EQ( col.total_length(), 0u );
    EXPECT_TRUE( col.fasta_files().empty() );
    EXPECT_TRUE( col.sequences().empty() );
    EXPECT_EQ( col.begin(), col.end() );
}

TEST(ReferenceCollection, AddSingleSequence)
{
    ReferenceCollection col;
    col.add( "ref.fa", "chr1 some description", 1000u );

    EXPECT_FALSE( col.empty() );
    EXPECT_EQ( col.size(), 1u );
    EXPECT_EQ( col.total_length(), 1000u );
    EXPECT_EQ( col.global_offset_of(0), 0u );
    EXPECT_EQ( col[0].length, 1000u );
    EXPECT_EQ( col[0].header, "chr1 some description" );
    EXPECT_EQ( col[0].file_index, 0u );
    ASSERT_EQ( col.fasta_files().size(), 1u );
    EXPECT_EQ( col.fasta_file(0), "ref.fa" );
}

TEST(ReferenceCollection, AddMultipleSequences)
{
    ReferenceCollection col;
    col.add( "ref.fa", "chr1", 100u );
    col.add( "ref.fa", "chr2", 200u );
    col.add( "ref.fa", "chr3",  50u );

    EXPECT_EQ( col.size(), 3u );
    EXPECT_EQ( col.total_length(), 350u );

    EXPECT_EQ( col.global_offset_of(0),   0u );
    EXPECT_EQ( col.global_offset_of(1), 100u );
    EXPECT_EQ( col.global_offset_of(2), 300u );

    EXPECT_EQ( col[0].length, 100u );
    EXPECT_EQ( col[1].length, 200u );
    EXPECT_EQ( col[2].length,  50u );
}

TEST(ReferenceCollection, MultipleFiles)
{
    ReferenceCollection col;
    col.add( "a.fa", "chr1", 10u );
    col.add( "b.fa", "chr2", 20u );
    col.add( "b.fa", "chr3", 30u );

    ASSERT_EQ( col.fasta_files().size(), 2u );
    EXPECT_EQ( col.fasta_file(0), "a.fa" );
    EXPECT_EQ( col.fasta_file(1), "b.fa" );

    EXPECT_EQ( col[0].file_index, 0u );
    EXPECT_EQ( col[1].file_index, 1u );
    EXPECT_EQ( col[2].file_index, 1u );
}

TEST(ReferenceCollection, SameFileTwiceIsDeduped)
{
    ReferenceCollection col;
    col.add( "same.fa", "seq1", 10u );
    col.add( "same.fa", "seq2", 20u );

    ASSERT_EQ( col.fasta_files().size(), 1u );
    EXPECT_EQ( col[0].file_index, 0u );
    EXPECT_EQ( col[1].file_index, 0u );
}

// =================================================================================================
//     SequenceRecord::name()
// =================================================================================================

TEST(ReferenceCollection, SequenceRecordName)
{
    ReferenceCollection col;
    col.add( "f.fa", "plain",               10u );
    col.add( "f.fa", "with space and more", 10u );
    col.add( "f.fa", "with\ttab",           10u );
    col.add( "f.fa", "multi word one two",  10u );

    EXPECT_EQ( col[0].name(), "plain" );
    EXPECT_EQ( col[1].name(), "with" );
    EXPECT_EQ( col[2].name(), "with" );
    EXPECT_EQ( col[3].name(), "multi" );
}

// =================================================================================================
//     add() overloads
// =================================================================================================

TEST(ReferenceCollection, SitesOverload)
{
    ReferenceCollection col;
    col.add( "f.fa", "chr1", std::string("ACGTACGT") );

    EXPECT_EQ( col.size(), 1u );
    EXPECT_EQ( col[0].length, 8u );
    EXPECT_EQ( col[0].sites, "ACGTACGT" );
    EXPECT_EQ( col.total_length(), 8u );
}

TEST(ReferenceCollection, MoveOverload)
{
    genesis::sequence::Sequence seq( "chr1 desc", "ACGTNNNN" );
    ReferenceCollection col;
    col.add( "f.fa", std::move(seq) );

    // sites were moved out of the original Sequence
    EXPECT_TRUE( seq.sites().empty() );
    EXPECT_EQ( col[0].length, 8u );
    EXPECT_EQ( col[0].sites, "ACGTNNNN" );
    EXPECT_EQ( col[0].header, "chr1 desc" );
}

TEST(ReferenceCollection, CopyOverload)
{
    genesis::sequence::Sequence seq( "chr1 desc", "ACGT" );
    ReferenceCollection col;
    col.add( "f.fa", seq );

    // original Sequence is untouched
    EXPECT_EQ( seq.sites(), "ACGT" );
    EXPECT_EQ( col[0].length, 4u );
    EXPECT_EQ( col[0].sites, "ACGT" );
}

TEST(ReferenceCollection, ZeroLengthThrows)
{
    ReferenceCollection col;
    EXPECT_THROW( col.add( "f.fa", "chr1", std::uint64_t(0) ), std::invalid_argument );
    EXPECT_THROW( col.add( "f.fa", "chr1", std::string("") ), std::invalid_argument );

    // Collection must remain empty after the failed adds.
    EXPECT_TRUE( col.empty() );
    EXPECT_EQ( col.total_length(), 0u );
}

// =================================================================================================
//     clear()
// =================================================================================================

TEST(ReferenceCollection, ClearAndReuse)
{
    ReferenceCollection col;
    col.add( "f.fa", "chr1", 100u );
    col.add( "f.fa", "chr2", 200u );

    col.clear();
    EXPECT_TRUE( col.empty() );
    EXPECT_EQ( col.size(), 0u );
    EXPECT_EQ( col.total_length(), 0u );
    EXPECT_TRUE( col.fasta_files().empty() );

    col.add( "g.fa", "chr3", 50u );
    EXPECT_EQ( col.size(), 1u );
    EXPECT_EQ( col.total_length(), 50u );
    EXPECT_EQ( col.fasta_file(0), "g.fa" );
}

// =================================================================================================
//     Iterators
// =================================================================================================

TEST(ReferenceCollection, Iterator)
{
    ReferenceCollection col;
    col.add( "f.fa", "chr1", 10u );
    col.add( "f.fa", "chr2", 20u );
    col.add( "f.fa", "chr3", 30u );

    std::vector<std::string> names;
    for( auto const& rec : col ) {
        names.push_back( rec.name() );
    }

    ASSERT_EQ( names.size(), 3u );
    EXPECT_EQ( names[0], "chr1" );
    EXPECT_EQ( names[1], "chr2" );
    EXPECT_EQ( names[2], "chr3" );
}

// =================================================================================================
//     find_by_global_position()
//
//     Five sequences with lengths 100, 1, 200, 50, 300 → total 651.
//     Offsets: chr1 [0,100)  chr2 [100,101)  chr3 [101,301)  chr4 [301,351)  chr5 [351,651)
// =================================================================================================

TEST(ReferenceCollection, FindByGlobalPosition)
{
    ReferenceCollection col;
    col.add( "f.fa", "chr1", 100u );   // [0,   100)
    col.add( "f.fa", "chr2",   1u );   // [100, 101)  — single-base sequence
    col.add( "f.fa", "chr3", 200u );   // [101, 301)
    col.add( "f.fa", "chr4",  50u );   // [301, 351)
    col.add( "f.fa", "chr5", 300u );   // [351, 651)

    ASSERT_EQ( col.total_length(), 651u );

    // First position of the entire collection.
    {
        auto r = col.find_by_global_position(0);
        EXPECT_EQ( r.index, 0u );
        EXPECT_EQ( r.local_offset, 0u );
        EXPECT_EQ( r.sequence->name(), "chr1" );
    }

    // Middle of chr1.
    {
        auto r = col.find_by_global_position(50);
        EXPECT_EQ( r.index, 0u );
        EXPECT_EQ( r.local_offset, 50u );
    }

    // Last position of chr1, one before the boundary.
    {
        auto r = col.find_by_global_position(99);
        EXPECT_EQ( r.index, 0u );
        EXPECT_EQ( r.local_offset, 99u );
    }

    // First (and only) position of chr2 — boundary crossing from chr1.
    {
        auto r = col.find_by_global_position(100);
        EXPECT_EQ( r.index, 1u );
        EXPECT_EQ( r.local_offset, 0u );
        EXPECT_EQ( r.sequence->name(), "chr2" );
    }

    // First position of chr3 — boundary crossing from the single-base chr2.
    {
        auto r = col.find_by_global_position(101);
        EXPECT_EQ( r.index, 2u );
        EXPECT_EQ( r.local_offset, 0u );
        EXPECT_EQ( r.sequence->name(), "chr3" );
    }

    // Middle of chr3.
    {
        auto r = col.find_by_global_position(200);
        EXPECT_EQ( r.index, 2u );
        EXPECT_EQ( r.local_offset, 99u );   // 200 - 101 = 99
    }

    // Last position of chr3.
    {
        auto r = col.find_by_global_position(300);
        EXPECT_EQ( r.index, 2u );
        EXPECT_EQ( r.local_offset, 199u );  // 300 - 101 = 199
    }

    // First position of chr4.
    {
        auto r = col.find_by_global_position(301);
        EXPECT_EQ( r.index, 3u );
        EXPECT_EQ( r.local_offset, 0u );
        EXPECT_EQ( r.sequence->name(), "chr4" );
    }

    // Last position of chr4.
    {
        auto r = col.find_by_global_position(350);
        EXPECT_EQ( r.index, 3u );
        EXPECT_EQ( r.local_offset, 49u );   // 350 - 301 = 49
    }

    // First position of chr5.
    {
        auto r = col.find_by_global_position(351);
        EXPECT_EQ( r.index, 4u );
        EXPECT_EQ( r.local_offset, 0u );
        EXPECT_EQ( r.sequence->name(), "chr5" );
    }

    // Last position of the entire collection.
    {
        auto r = col.find_by_global_position(650);
        EXPECT_EQ( r.index, 4u );
        EXPECT_EQ( r.local_offset, 299u );  // 650 - 351 = 299
        EXPECT_EQ( r.sequence->name(), "chr5" );
    }
}

TEST(ReferenceCollection, FindByGlobalPositionErrors)
{
    // Empty collection throws on any position.
    {
        ReferenceCollection col;
        EXPECT_THROW( col.find_by_global_position(0), std::out_of_range );
    }

    // Position exactly at total_length() is out of range (half-open interval).
    {
        ReferenceCollection col;
        col.add( "f.fa", "chr1", 100u );
        EXPECT_NO_THROW( col.find_by_global_position(99) );
        EXPECT_THROW( col.find_by_global_position(100), std::out_of_range );
        EXPECT_THROW( col.find_by_global_position(101), std::out_of_range );
        EXPECT_THROW(
            col.find_by_global_position( std::numeric_limits<std::uint64_t>::max() ),
            std::out_of_range
        );
    }
}

TEST(ReferenceCollection, FindByGlobalPositionSingleBase)
{
    // Degenerate case: one sequence of length 1.
    ReferenceCollection col;
    col.add( "f.fa", "chr1", 1u );

    auto r = col.find_by_global_position(0);
    EXPECT_EQ( r.index, 0u );
    EXPECT_EQ( r.local_offset, 0u );
    EXPECT_THROW( col.find_by_global_position(1), std::out_of_range );
}

// =================================================================================================
//     JSON round-trip
// =================================================================================================

TEST(ReferenceCollectionJson, RoundTrip)
{
    ReferenceCollection original;
    original.add( "a.fa", "chr1 description", 1000u );
    original.add( "a.fa", "chr2",              500u );
    original.add( "b.fa", "chrX",             2000u );

    auto const clone = from_json( to_json(original) );

    ASSERT_EQ( clone.size(), 3u );
    ASSERT_EQ( clone.fasta_files().size(), 2u );

    EXPECT_EQ( clone.fasta_file(0), "a.fa" );
    EXPECT_EQ( clone.fasta_file(1), "b.fa" );

    EXPECT_EQ( clone[0].header,     "chr1 description" );
    EXPECT_EQ( clone[0].length,     1000u );
    EXPECT_EQ( clone[0].file_index, 0u );

    EXPECT_EQ( clone[1].header,     "chr2" );
    EXPECT_EQ( clone[1].length,     500u );
    EXPECT_EQ( clone[1].file_index, 0u );

    EXPECT_EQ( clone[2].header,     "chrX" );
    EXPECT_EQ( clone[2].length,     2000u );
    EXPECT_EQ( clone[2].file_index, 1u );

    EXPECT_EQ( clone.total_length(),       original.total_length() );
    EXPECT_EQ( clone.global_offset_of(0),  original.global_offset_of(0) );
    EXPECT_EQ( clone.global_offset_of(1),  original.global_offset_of(1) );
    EXPECT_EQ( clone.global_offset_of(2),  original.global_offset_of(2) );
}

TEST(ReferenceCollectionJson, SitesNotStored)
{
    ReferenceCollection original;
    original.add( "f.fa", "chr1", std::string("ACGTACGT") );

    auto const clone = from_json( to_json(original) );
    EXPECT_EQ( clone[0].length, 8u );
    EXPECT_TRUE( clone[0].sites.empty() );
}

TEST(ReferenceCollectionJson, EmptyCollection)
{
    ReferenceCollection original;
    auto const clone = from_json( to_json(original) );
    EXPECT_TRUE( clone.empty() );
    EXPECT_EQ( clone.total_length(), 0u );
}

// =================================================================================================
//     build_collection()
// =================================================================================================

TEST(BuildCollection, FallbackScan)
{
    TempFile fasta;
    {
        std::ofstream f( fasta.path );
        f << ">chr1 some description\nACGTACGT\n";
        f << ">chr2\nGGGG\n";
    }

    auto const col = build_collection({ fasta.path });

    ASSERT_EQ( col.size(), 2u );
    EXPECT_EQ( col[0].name(), "chr1" );
    EXPECT_EQ( col[0].length, 8u );
    EXPECT_EQ( col[1].name(), "chr2" );
    EXPECT_EQ( col[1].length, 4u );
    EXPECT_EQ( col.total_length(), 12u );
    EXPECT_EQ( col.fasta_file(0), fasta.path );
}

TEST(BuildCollection, FaiPath)
{
    TempFile fasta;
    {
        std::ofstream f( fasta.path );
        // Minimal FASTA content; build_collection won't read this when the FAI is found.
        f << ">placeholder\nAAAA\n";
    }

    // Write a .fai with lengths that intentionally differ from the FASTA above,
    // proving that build_collection reads lengths from the FAI and ignores the FASTA.
    std::string const fai_path = fasta.path + ".fai";
    {
        std::ofstream fai( fai_path );
        // Columns: name, length, offset, line_bases, line_width (last three are unused by our code)
        fai << "chr1\t1000\t6\t1000\t1001\n";
        fai << "chr2\t2000\t1013\t2000\t2001\n";
    }

    auto const col = build_collection({ fasta.path });
    ::unlink( fai_path.c_str() );

    ASSERT_EQ( col.size(), 2u );
    EXPECT_EQ( col[0].name(), "chr1" );
    EXPECT_EQ( col[0].length, 1000u );
    EXPECT_EQ( col[1].name(), "chr2" );
    EXPECT_EQ( col[1].length, 2000u );
    EXPECT_EQ( col.total_length(), 3000u );
}

TEST(BuildCollection, StoreSites)
{
    TempFile fasta;
    {
        std::ofstream f( fasta.path );
        f << ">chr1\nACGT\n>chr2\nGGCC\n";
    }

    auto const col = build_collection({ fasta.path }, /*store_sites=*/true);

    ASSERT_EQ( col.size(), 2u );
    EXPECT_EQ( col[0].sites, "ACGT" );
    EXPECT_EQ( col[1].sites, "GGCC" );
}
