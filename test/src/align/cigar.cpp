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
 * @brief Tests for CIGAR utility functions (cigar.hpp).
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/align/cigar.hpp"

#include <stdexcept>
#include <vector>

using namespace spear::align;

// BAM uint32_t CIGAR encoding helpers: (length << 4) | op_code
// Op codes: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
static constexpr uint32_t op( uint32_t len, uint32_t code ) { return (len << 4) | code; }

// =================================================================================================
//     cigar_to_string
// =================================================================================================

TEST( CigarToString, Empty )
{
    EXPECT_EQ( cigar_to_string( {} ), "" );
}

TEST( CigarToString, SingleMatch )
{
    EXPECT_EQ( cigar_to_string( { op(8, 7) } ), "8=" );
}

TEST( CigarToString, AllOps )
{
    // M I D N S H P = X
    std::vector<uint32_t> cig = {
        op(1,0), op(2,1), op(3,2), op(4,3),
        op(5,4), op(6,5), op(7,6), op(8,7), op(9,8)
    };
    EXPECT_EQ( cigar_to_string( cig ), "1M2I3D4N5S6H7P8=9X" );
}

TEST( CigarToString, TypicalAlignment )
{
    // 30= 2X 5= 1I 10=
    std::vector<uint32_t> cig = { op(30,7), op(2,8), op(5,7), op(1,1), op(10,7) };
    EXPECT_EQ( cigar_to_string( cig ), "30=2X5=1I10=" );
}

// =================================================================================================
//     edit_distance_from_cigar (single-argument, X/= CIGAR only)
// =================================================================================================

TEST( EditDistanceSingleArg, PerfectMatch )
{
    // All = ops: edit distance is 0
    EXPECT_EQ( edit_distance_from_cigar( { op(100, 7) } ), 0 );
}

TEST( EditDistanceSingleArg, Mismatches )
{
    // 8= 2X 5= → 2 mismatches
    std::vector<uint32_t> cig = { op(8,7), op(2,8), op(5,7) };
    EXPECT_EQ( edit_distance_from_cigar( cig ), 2 );
}

TEST( EditDistanceSingleArg, Insertion )
{
    // 10= 3I 5= → 3 inserted bases
    std::vector<uint32_t> cig = { op(10,7), op(3,1), op(5,7) };
    EXPECT_EQ( edit_distance_from_cigar( cig ), 3 );
}

TEST( EditDistanceSingleArg, Deletion )
{
    // 10= 2D 5= → 2 deleted bases
    std::vector<uint32_t> cig = { op(10,7), op(2,2), op(5,7) };
    EXPECT_EQ( edit_distance_from_cigar( cig ), 2 );
}

TEST( EditDistanceSingleArg, Mixed )
{
    // 5= 1X 3= 2I 4= 1D 6= → NM = 1+2+1 = 4
    std::vector<uint32_t> cig = { op(5,7), op(1,8), op(3,7), op(2,1), op(4,7), op(1,2), op(6,7) };
    EXPECT_EQ( edit_distance_from_cigar( cig ), 4 );
}

TEST( EditDistanceSingleArg, ThrowsOnM )
{
    std::vector<uint32_t> cig = { op(10, 0) };
    EXPECT_THROW( edit_distance_from_cigar( cig ), std::invalid_argument );
}

TEST( EditDistanceSingleArg, ThrowsOnUnexpectedOp )
{
    // N(3), S(4), H(5), P(6) are all unexpected from WFA2 gap-affine
    for( uint32_t bad : { 3u, 4u, 5u, 6u } ) {
        EXPECT_THROW( edit_distance_from_cigar( { op(1, bad) } ), std::invalid_argument )
            << "op code " << bad << " should throw";
    }
}

// =================================================================================================
//     edit_distance_from_cigar (two-sequence overload, handles both M and X/= CIGAR)
// =================================================================================================

TEST( EditDistanceTwoSeq, XeqCigar )
{
    // With X/= CIGAR the sequences are present but not consulted for M ops; result must
    // match the single-arg overload exactly.
    //   query:  ACGTACGT
    //   ref:    ACGTACGT  (passed as window; ref_begin=0)
    //   CIGAR:  6= 1X 1=   → NM=1
    std::string query   = "ACGTACXT";  // position 6 differs
    std::string ref_win = "ACGTACGT";
    std::vector<uint32_t> cig = { op(6,7), op(1,8), op(1,7) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 0 ), 1 );
    EXPECT_EQ( edit_distance_from_cigar( cig ), 1 );  // single-arg agrees
}

TEST( EditDistanceTwoSeq, McigarPerfectMatch )
{
    // M CIGAR, all bases match → NM=0
    std::string query   = "ACGT";
    std::string ref_win = "ACGT";
    std::vector<uint32_t> cig = { op(4, 0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 0 ), 0 );
}

TEST( EditDistanceTwoSeq, McigarMismatches )
{
    // M CIGAR: query ACGT vs ref ACCT → 1 mismatch at position 2
    std::string query   = "ACGT";
    std::string ref_win = "ACCT";
    std::vector<uint32_t> cig = { op(4, 0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 0 ), 1 );
}

TEST( EditDistanceTwoSeq, McigarWithIndels )
{
    // M CIGAR: 3M 1I 3M against ref that matches the M spans
    // query:  A C G A C G T  (length 7: ACG + I:A + CGT but I is query-only)
    // Actually let's construct this carefully:
    // query:  "ACGACGT" (7 bases), CIGAR: 3M 1I 3M
    //   M span 1: query[0..2]="ACG" vs ref[0..2]="ACG" → 0 mismatches
    //   I span:   query[3]="A" consumed, no ref advance
    //   M span 2: query[4..6]="CGT" vs ref[3..5]="CGT" → 0 mismatches
    // NM = 0 matches + 1 insertion = 1
    std::string query   = "ACGACGT";
    std::string ref_win = "ACGCGT";
    std::vector<uint32_t> cig = { op(3,0), op(1,1), op(3,0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 0 ), 1 );
}

TEST( EditDistanceTwoSeq, McigarWithDeletion )
{
    // query: "ACGCGT" (6 bases), CIGAR: 3M 1D 3M
    //   M span 1: query[0..2]="ACG" vs ref[0..2]="ACG" → 0 mismatches
    //   D span:   ref[3]="A" consumed, no query advance → 1 deletion
    //   M span 2: query[3..5]="CGT" vs ref[4..6]="CGT" → 0 mismatches
    // NM = 0 mismatches + 1 deletion = 1
    std::string query   = "ACGCGT";
    std::string ref_win = "ACGACGT";
    std::vector<uint32_t> cig = { op(3,0), op(1,2), op(3,0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 0 ), 1 );
}

TEST( EditDistanceTwoSeq, RefBeginOffset )
{
    // Alignment starts at ref_begin=4 inside a longer window.
    // query: "ACGT", ref_window: "NNNNACGTNNN", ref_begin=4
    // CIGAR: 4M, all bases match → NM=0
    std::string query   = "ACGT";
    std::string ref_win = "NNNNACGTNNN";
    std::vector<uint32_t> cig = { op(4, 0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 4 ), 0 );
}

TEST( EditDistanceTwoSeq, RefBeginOffsetWithMismatch )
{
    // Same setup but position 2 in query differs from ref_window[4+2]=ref_window[6]
    // query: "ACXT", ref_window: "NNNNACGTNNN" → position 2: X vs G → 1 mismatch
    std::string query   = "ACXT";
    std::string ref_win = "NNNNACGTNNN";
    std::vector<uint32_t> cig = { op(4, 0) };
    EXPECT_EQ( edit_distance_from_cigar( cig, query, ref_win, 4 ), 1 );
}

TEST( EditDistanceTwoSeq, ThrowsOnUnexpectedOp )
{
    std::string q = "ACGT", r = "ACGT";
    for( uint32_t bad : { 3u, 4u, 5u, 6u } ) {
        EXPECT_THROW( edit_distance_from_cigar( { op(1, bad) }, q, r, 0 ), std::invalid_argument )
            << "op code " << bad << " should throw";
    }
}
