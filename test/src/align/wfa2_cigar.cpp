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
 * @brief Unit tests for build_sam_cigar_from_wfa2_ops() (wfa2_cigar.hpp), exercised
 * directly with synthetic WFA2 op strings -- no real WFA2 alignment is run here.
 *
 * See test/src/align/wfa2.cpp for integration-level tests that go through the real
 * Wfa2Aligner.
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/align/cigar.hpp"
#include "spear/align/wfa2_cigar.hpp"

using namespace spear::align;

// =================================================================================================
//     Basic matches
// =================================================================================================

TEST( Wfa2CigarBuild, PerfectMatchExtended )
{
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMM", "ACGT", "ACGT", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2CigarBuild, PerfectMatchCollapsed )
{
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMM", "ACGT", "ACGT", false );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4M" );
}

TEST( Wfa2CigarBuild, Empty )
{
    auto const r = build_sam_cigar_from_wfa2_ops( "", "", "", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   0 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_TRUE( r.cigar.empty() );
}

// =================================================================================================
//     Real (WFA2-charged) mismatches: raw 'X'
// =================================================================================================

TEST( Wfa2CigarBuild, RealMismatchExtended )
{
    // query: ACGT vs target: ACCT -- position 2 is a real (scored) mismatch.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMXM", "ACGT", "ACCT", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "2=1X1=" );
}

TEST( Wfa2CigarBuild, RealMismatchCollapsed )
{
    auto const r = build_sam_cigar_from_wfa2_ops( "MMXM", "ACGT", "ACCT", false );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4M" );  // X collapses into the surrounding M runs
}

// =================================================================================================
//     'M' with mismatched characters: the aDNA damage-tolerance case
// =================================================================================================
//
// WFA2 only ever emits 'M' for a position where characters differ if the match function
// used during alignment (e.g. the damage-aware one in wfa2.cpp) accepted it as a tolerated
// substitution. Such positions must still be reported as mismatches in the CIGAR,
// even though WFA2's own raw op label says "match".

TEST( Wfa2CigarBuild, DamageTolerantMismatchExtended )
{
    // query: ACTGT vs target: ACCGT -- position 2 (T vs C) is WFA2-labelled 'M' (accepted
    // as a match during alignment) but the characters actually differ.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMMM", "ACTGT", "ACCGT", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   5 );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "2=1X2=" );  // NOT "5=" -- the actual bug
}

TEST( Wfa2CigarBuild, DamageTolerantMismatchCollapsed )
{
    // Same scenario, collapsed CIGAR: the emitted op is still plain M either way, but
    // edit_distance must still correctly detect the mismatch.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMMM", "ACTGT", "ACCGT", false );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "5M" );
}

TEST( Wfa2CigarBuild, RealisticCT5PrimeDamage )
{
    // Mirrors an actual C->T 5' deamination: read="TACGT" vs ref="CACGT", position 0.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMMM", "TACGT", "CACGT", true );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "1X4=" );
}

// =================================================================================================
//     Leading / trailing reference padding (free-end clips)
// =================================================================================================

TEST( Wfa2CigarBuild, LeadingClip )
{
    // 2 bases of reference padding before the read starts aligning.
    auto const r = build_sam_cigar_from_wfa2_ops( "IIMMMM", "ACGT", "XXACGT", true );
    EXPECT_EQ( r.ref_begin, 2 );
    EXPECT_EQ( r.ref_end,   6 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2CigarBuild, TrailingClip )
{
    // 2 bases of reference padding after the read finishes aligning.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMMII", "ACGT", "ACGTXX", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2CigarBuild, BothClips )
{
    auto const r = build_sam_cigar_from_wfa2_ops( "IMMMMI", "ACGT", "XACGTX", true );
    EXPECT_EQ( r.ref_begin, 1 );
    EXPECT_EQ( r.ref_end,   5 );
    EXPECT_EQ( r.edit_distance, 0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

// =================================================================================================
//     Indels: WFA2's I/D convention is swapped relative to SAM
// =================================================================================================

TEST( Wfa2CigarBuild, InsertionFromRawD )
{
    // Raw 'D' (query-consuming only) must become SAM 'I'.
    // query: ACAGT (extra 'A' at position 2) vs target: ACGT.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMDMM", "ACAGT", "ACGT", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "2=1I2=" );
}

TEST( Wfa2CigarBuild, DeletionFromRawI )
{
    // Raw 'I' (reference-consuming only) in the middle of the alignment -- not at either
    // boundary, so it's a real deletion, not free-end padding -- must become SAM 'D'.
    // query: ACGT vs target: ACXGT (extra 'X' in the reference at position 2).
    auto const r = build_sam_cigar_from_wfa2_ops( "MMIMM", "ACGT", "ACXGT", true );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   5 );
    EXPECT_EQ( r.edit_distance, 1 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "2=1D2=" );
}

TEST( Wfa2CigarBuild, IndelsUnaffectedByExtendedFlag )
{
    // I/D are unambiguous regardless of use_extended_cigar -- only M/X collapse.
    auto const extended  = build_sam_cigar_from_wfa2_ops( "MMDMM", "ACAGT", "ACGT", true );
    auto const collapsed = build_sam_cigar_from_wfa2_ops( "MMDMM", "ACAGT", "ACGT", false );
    EXPECT_EQ( cigar_to_string( extended.cigar  ), "2=1I2=" );
    EXPECT_EQ( cigar_to_string( collapsed.cigar ), "2M1I2M" );
    EXPECT_EQ( extended.edit_distance, collapsed.edit_distance );
}

// =================================================================================================
//     Run merging
// =================================================================================================

TEST( Wfa2CigarBuild, AdjacentSameOpsMerge )
{
    // Two separate mismatches with a match run between them must still merge correctly
    // with their neighbours and not merge across the mismatch.
    // query: AXCXA vs target: AACAA -- positions 1 and 3 differ.
    auto const r = build_sam_cigar_from_wfa2_ops( "MMMMM", "AXCXA", "AACAA", true );
    EXPECT_EQ( r.edit_distance, 2 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "1=1X1=1X1=" );
}
