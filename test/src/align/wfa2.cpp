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
 * @brief Tests for Wfa2Aligner and surrounding utilities (wfa2.hpp, cigar.hpp).
 *
 * @file
 * @ingroup test
 */

#include "src/common.hpp"

#include "spear/align/cigar.hpp"
#include "spear/align/wfa2.hpp"

#include <stdexcept>
#include <string>

using namespace spear::align;

// Default penalty settings used throughout: mismatch=4, gap_open=4, gap_extend=2.
// Score convention: 0 = perfect, more negative = worse.
// A single mismatch costs -4; a 1-base gap costs -(gap_open + gap_extend) = -6.

// =================================================================================================
//     Settings validation
// =================================================================================================

TEST( Wfa2Settings, NegativeMismatchThrows )
{
    Wfa2Settings s;
    s.mismatch = -1;
    EXPECT_THROW( Wfa2Aligner{ s }, std::invalid_argument );
}

TEST( Wfa2Settings, NegativeGapOpenThrows )
{
    Wfa2Settings s;
    s.gap_open = -1;
    EXPECT_THROW( Wfa2Aligner{ s }, std::invalid_argument );
}

TEST( Wfa2Settings, NegativeGapExtendThrows )
{
    Wfa2Settings s;
    s.gap_extend = -1;
    EXPECT_THROW( Wfa2Aligner{ s }, std::invalid_argument );
}

TEST( Wfa2Settings, NegativeMaxStepsThrows )
{
    Wfa2Settings s;
    s.max_steps = -1;
    EXPECT_THROW( Wfa2Aligner{ s }, std::invalid_argument );
}

// =================================================================================================
//     Wrong-scope guards
// =================================================================================================

TEST( Wfa2Scope, AlignScoreOnCigarAlignerThrows )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };
    EXPECT_THROW( aligner.align_score( "ACGT", "ACGT" ), std::logic_error );
}

TEST( Wfa2Scope, AlignCigarOnScoreAlignerThrows )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };
    EXPECT_THROW( aligner.align_cigar( "ACGT", "ACGT" ), std::logic_error );
}

// =================================================================================================
//     Empty / degenerate inputs
// =================================================================================================

TEST( Wfa2Empty, EmptyQueryScore )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "", "ACGT" );
    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   0 );
}

TEST( Wfa2Empty, EmptyQueryCigar )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "", "ACGT" );
    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   0 );
    EXPECT_TRUE( r.cigar.empty() );
}

TEST( Wfa2Empty, BothEmptyScore )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "", "" );
    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   0 );
}

TEST( Wfa2Empty, BothEmptyCigar )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "", "" );
    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   0 );
    EXPECT_TRUE( r.cigar.empty() );
}

TEST( Wfa2Empty, EmptyRefScoreThrows )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };
    EXPECT_THROW( aligner.align_score( "ACGT", "" ), std::invalid_argument );
}

TEST( Wfa2Empty, EmptyRefCigarThrows )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };
    EXPECT_THROW( aligner.align_cigar( "ACGT", "" ), std::invalid_argument );
}

// =================================================================================================
//     align_score
// =================================================================================================

TEST( Wfa2AlignScore, PerfectMatch )
{
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "ACGTACGT", "ACGTACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
    EXPECT_EQ( r.ref_end, 8 );
}

TEST( Wfa2AlignScore, SingleMismatch )
{
    // Position 4: A vs T → penalty -4
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "ACGTTCGT", "ACGTACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  -4 );
    EXPECT_EQ( r.ref_end, 8 );
}

TEST( Wfa2AlignScore, FittingInMiddle )
{
    // Read aligns perfectly to the middle of the reference window; ref_end should be 8 (not 12).
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "ACGT", "TTTTACGTTTTT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
    EXPECT_EQ( r.ref_end, 8 );
}

TEST( Wfa2AlignScore, ScoreCap )
{
    // Very low max_steps forces early termination on a highly divergent alignment.
    Wfa2Settings s;
    s.scope     = Wfa2Settings::Scope::Score;
    s.max_steps = 1;
    Wfa2Aligner aligner{ s };

    // All mismatches: AAAA vs TTTT should exceed max_steps=1
    auto const r = aligner.align_score( "AAAA", "TTTT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::ScoreCapExceeded );
}

// =================================================================================================
//     align_cigar
// =================================================================================================

TEST( Wfa2AlignCigar, PerfectMatchExtendedCigar )
{
    // Default use_extended_cigar=true → CIGAR should be all = (op 7)
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    std::string const query = "ACGTACGT";
    auto const r = aligner.align_cigar( query, query );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "8=" );
}

TEST( Wfa2AlignCigar, MCigarFormat )
{
    // use_extended_cigar=false → same alignment should produce M ops instead of =
    Wfa2Settings s;
    s.scope              = Wfa2Settings::Scope::Cigar;
    s.use_extended_cigar = false;
    Wfa2Aligner aligner{ s };

    std::string const query = "ACGTACGT";
    auto const r = aligner.align_cigar( query, query );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "8M" );
}

TEST( Wfa2AlignCigar, SingleMismatch )
{
    // query: ACGTTCGT, ref: ACGTACGT → position 4: T vs A → X in extended CIGAR
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGTTCGT", "ACGTACGT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     -4 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=1X3=" );
}

TEST( Wfa2AlignCigar, Insertion )
{
    // query: ACGAACGT (8 bases), ref: ACGACGT (7 bases)
    // Expected: 3= 1I 4= score=-(gap_open+gap_extend)=-6
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGAACGT", "ACGACGT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     -6 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   7 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=1I3=" );
}

TEST( Wfa2AlignCigar, Deletion )
{
    // query: ACGACGT (7 bases), ref: ACGAACGT (8 bases)
    // Expected: 4= 1D 3= score=-6 (WFA2 right-aligns indels)
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGACGT", "ACGAACGT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     -6 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=1D3=" );
}

TEST( Wfa2AlignCigar, FittingInMiddle )
{
    // Read aligns to positions 4–7 of the reference window; flanks don't match.
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGT", "TTTTACGTTTTT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 4 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2AlignCigar, FittingClipEnd )
{
    // Read aligns flush to the start of the window; clip only at the trailing end.
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGT", "ACGTTTTT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 0 );
    EXPECT_EQ( r.ref_end,   4 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2AlignCigar, FittingClipStart )
{
    // Read aligns flush to the end of the window; clip only at the leading end.
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGT", "TTTTACGT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     0 );
    EXPECT_EQ( r.ref_begin, 4 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "4=" );
}

TEST( Wfa2AlignCigar, FittingBothClipsWithMismatch )
{
    // Read aligns in the middle of the window and has a mismatch in the core alignment.
    // query: ACXT, ref window: TTTTACGTTTTT → core ref region [4,8) = ACGT, pos 2 differs.
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACXT", "TTTTACGTTTTT" );

    EXPECT_EQ( r.status,    Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,     -4 );
    EXPECT_EQ( r.ref_begin, 4 );
    EXPECT_EQ( r.ref_end,   8 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "2=1X1=" );
}

TEST( Wfa2AlignCigar, ScoreCap )
{
    Wfa2Settings s;
    s.scope     = Wfa2Settings::Scope::Cigar;
    s.max_steps = 1;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "AAAA", "TTTT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::ScoreCapExceeded );
    EXPECT_TRUE( r.cigar.empty() );
}

TEST( Wfa2AlignCigar, ScoreAgreesBetweenModes )
{
    // align_score and align_cigar must report the same score and ref_end for identical inputs.
    std::string const query  = "ACGTTCGT";
    std::string const target = "ACGTACGT";

    Wfa2Settings ss;
    ss.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner score_aligner{ ss };

    Wfa2Settings cs;
    cs.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner cigar_aligner{ cs };

    auto const rs = score_aligner.align_score( query, target );
    auto const rc = cigar_aligner.align_cigar( query, target );

    EXPECT_EQ( rs.status,  Wfa2Result::Status::Ok );
    EXPECT_EQ( rc.status,  Wfa2Result::Status::Ok );
    EXPECT_EQ( rs.score,   rc.score );
    EXPECT_EQ( rs.ref_end, rc.ref_end );
}

// =================================================================================================
//     aDNA damage tolerance
// =================================================================================================

TEST( Wfa2Damage, CT5PrimeTreatedAsMatch )
{
    // Position 0: read=T, ref=C → C→T deamination damage; should score as perfect match.
    Wfa2Settings s;
    s.scope         = Wfa2Settings::Scope::Score;
    s.damage_ct_end = 3;  // zone covers positions 0, 1, 2
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "TACGT", "CACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
}

TEST( Wfa2Damage, CT5PrimeOutsideZonePenalised )
{
    // Same C→T substitution but at position 3, outside the ct_end=2 zone → regular mismatch.
    Wfa2Settings s;
    s.scope         = Wfa2Settings::Scope::Score;
    s.damage_ct_end = 2;  // zone covers positions 0, 1 only
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "ACGTA", "ACGCA" );  // pos 3: T vs C, outside zone
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  -4 );
}

TEST( Wfa2Damage, GA3PrimeTreatedAsMatch )
{
    // Positions 3 and 4 (last two of a 5-base read): read=A, ref=G → G→A damage.
    Wfa2Settings s;
    s.scope          = Wfa2Settings::Scope::Score;
    s.damage_ga_start = 3;  // zone covers positions 3 and 4
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "ACGAA", "ACGGG" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
}

TEST( Wfa2Damage, GA3PrimeOutsideZonePenalised )
{
    // G→A at position 1, but ga_start=3 means only positions 3+ are covered.
    Wfa2Settings s;
    s.scope           = Wfa2Settings::Scope::Score;
    s.damage_ga_start = 3;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "AAGGG", "AGGGG" );  // pos 1: A vs G, outside zone
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  -4 );
}

TEST( Wfa2Damage, NoDamageWithoutSettings )
{
    // Same C→T substitution at position 0, but no damage settings → penalised normally.
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Score;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_score( "TACGT", "CACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  -4 );
}

// The Score-only tests above only check that damage tolerance affects the alignment
// *score*. The tests below test that a tolerated substitution must still score as a match
// (score == 0), but must be reported as a mismatch (X, not =) in the CIGAR, and must be
// counted in edit_distance.

TEST( Wfa2Damage, CT5PrimeCigarShowsMismatch )
{
    Wfa2Settings s;
    s.scope         = Wfa2Settings::Scope::Cigar;
    s.damage_ct_end = 3;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "TACGT", "CACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );  // tolerated during alignment/scoring
    EXPECT_EQ( cigar_to_string( r.cigar ), "1X4=" );  // but reported as a mismatch in the CIGAR
    EXPECT_EQ( r.edit_distance, 1 );
}

TEST( Wfa2Damage, GA3PrimeCigarShowsMismatch )
{
    Wfa2Settings s;
    s.scope           = Wfa2Settings::Scope::Cigar;
    s.damage_ga_start = 3;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGAA", "ACGGG" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "3=2X" );
    EXPECT_EQ( r.edit_distance, 2 );
}

TEST( Wfa2Damage, CigarShowsMismatchWithCollapsedMCigar )
{
    // Same as CT5PrimeCigarShowsMismatch, but with use_extended_cigar=false: the CIGAR
    // itself is ambiguous (all M), but edit_distance must still detect the substitution.
    Wfa2Settings s;
    s.scope              = Wfa2Settings::Scope::Cigar;
    s.damage_ct_end      = 3;
    s.use_extended_cigar = false;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "TACGT", "CACGT" );
    EXPECT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( r.score,  0 );
    EXPECT_EQ( cigar_to_string( r.cigar ), "5M" );
    EXPECT_EQ( r.edit_distance, 1 );
}

// =================================================================================================
//     Edit distance cross-validation
// =================================================================================================

TEST( Wfa2EditDistance, XeqCigarMatchesExpectedNM )
{
    // Align with known mismatches and indels; check edit_distance_from_cigar agrees.
    // query: ACGTTCGT (1 mismatch at pos 4 vs ACGTACGT) → NM should be 1
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( "ACGTTCGT", "ACGTACGT" );
    ASSERT_EQ( r.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( edit_distance_from_cigar( r.cigar ), 1 );
}

TEST( Wfa2EditDistance, MCigarMatchesTwoSeqComputation )
{
    // With M CIGAR the single-arg overload throws; the two-sequence overload must agree
    // with the known expected NM.
    std::string const query  = "ACGTTCGT";
    std::string const target = "ACGTACGT";

    Wfa2Settings s;
    s.scope              = Wfa2Settings::Scope::Cigar;
    s.use_extended_cigar = false;
    Wfa2Aligner aligner{ s };

    auto const r = aligner.align_cigar( query, target );
    ASSERT_EQ( r.status, Wfa2Result::Status::Ok );

    // Single-arg must throw for M CIGAR
    EXPECT_THROW( edit_distance_from_cigar( r.cigar ), std::invalid_argument );

    // Two-sequence overload resolves the M ops and counts the one mismatch
    EXPECT_EQ( edit_distance_from_cigar( r.cigar, query, target, r.ref_begin ), 1 );
}

TEST( Wfa2EditDistance, InsertionAndDeletion )
{
    // Insertion test: NM should be 1 (the inserted base)
    Wfa2Settings s;
    s.scope = Wfa2Settings::Scope::Cigar;
    Wfa2Aligner aligner{ s };

    auto const ri = aligner.align_cigar( "ACGAACGT", "ACGACGT" );
    ASSERT_EQ( ri.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( edit_distance_from_cigar( ri.cigar ), 1 );

    // Deletion test: NM should also be 1 (the deleted base)
    auto const rd = aligner.align_cigar( "ACGACGT", "ACGAACGT" );
    ASSERT_EQ( rd.status, Wfa2Result::Status::Ok );
    EXPECT_EQ( edit_distance_from_cigar( rd.cigar ), 1 );
}
