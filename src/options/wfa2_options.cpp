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

#include "options/wfa2_options.hpp"

#include "tools/misc.hpp"

// =================================================================================================
//     Setup
// =================================================================================================

void Wfa2Options::add_wfa2_opts_to_app(
    CLI::App* sub, std::string const& group
) {
    internal_check(
        !mismatch_penalty.option && !gap_open_penalty.option && !gap_extend_penalty.option &&
        !damage_ct_5p_reach.option && !damage_ga_3p_reach.option && !damage_reach.option &&
        !max_steps.option && !basic_cigar.option,
        "Cannot use the same Wfa2Options object multiple times."
    );

    // Mismatch penalty.
    mismatch_penalty = sub->add_option(
        "--align-mismatch-penalty",
        mismatch_penalty.value,
        "Penalty for a mismatched base during alignment."
    );
    mismatch_penalty.option->check( CLI::PositiveNumber );
    mismatch_penalty.option->group( group );

    // Gap-open penalty.
    gap_open_penalty = sub->add_option(
        "--align-gap-open-penalty",
        gap_open_penalty.value,
        "Penalty for opening a gap (insertion or deletion) during alignment, in addition to "
        "--align-gap-extend-penalty for the first base of the gap."
    );
    gap_open_penalty.option->check( CLI::NonNegativeNumber );
    gap_open_penalty.option->group( group );

    // Gap-extension penalty.
    gap_extend_penalty = sub->add_option(
        "--align-gap-extend-penalty",
        gap_extend_penalty.value,
        "Penalty for each base of a gap (insertion or deletion) during alignment, following the "
        "gap-affine Smith-Waterman-Gotoh cost model."
    );
    gap_extend_penalty.option->check( CLI::PositiveNumber );
    gap_extend_penalty.option->group( group );

    // aDNA 5' damage tolerance.
    damage_ct_5p_reach = sub->add_option(
        "--align-damage-ct-5p-reach",
        damage_ct_5p_reach.value,
        "Number of bases from the 5' end of each read within which a C-to-T substitution "
        "(relative to the reference) is treated as a free match instead of a mismatch, to "
        "tolerate the characteristic aDNA cytosine deamination pattern. The default of 0 "
        "disables 5' damage tolerance."
    );
    damage_ct_5p_reach.option->check( CLI::NonNegativeNumber );
    damage_ct_5p_reach.option->group( group );

    // aDNA 3' damage tolerance.
    damage_ga_3p_reach = sub->add_option(
        "--align-damage-ga-3p-reach",
        damage_ga_3p_reach.value,
        "Number of bases from the 3' end of each read within which a G-to-A substitution "
        "(relative to the reference) is treated as a free match instead of a mismatch, to "
        "tolerate the characteristic aDNA cytosine deamination pattern. The default of 0 "
        "disables 3' damage tolerance."
    );
    damage_ga_3p_reach.option->check( CLI::NonNegativeNumber );
    damage_ga_3p_reach.option->group( group );

    // aDNA damage tolerance convenience option: sets both zones at once.
    damage_reach = sub->add_option(
        "--align-damage-reach",
        damage_reach.value,
        "Convenience option that sets --align-damage-ct-5p-reach and --align-damage-ga-3p-reach "
        "to the same value at once, for the common case of wanting the same aDNA damage tolerance "
        "at both ends of the read."
    );
    damage_reach.option->check( CLI::NonNegativeNumber );
    damage_reach.option->group( group );
    damage_reach.option->excludes( damage_ct_5p_reach.option );
    damage_reach.option->excludes( damage_ga_3p_reach.option );

    // Max alignment steps.
    max_steps = sub->add_option(
        "--align-max-steps",
        max_steps.value,
        "Safety cap on the alignment score (in WFA steps) before aborting a pathological "
        "alignment. The default is generous enough to guarantee completion at any realistic "
        "divergence for typical short-read lengths."
    );
    max_steps.option->check( CLI::NonNegativeNumber );
    max_steps.option->group( group );

    // Basic (M-only) CIGAR.
    basic_cigar = sub->add_flag(
        "--align-basic-cigar",
        basic_cigar.value,
        "If set, collapse the '=' (match) and 'X' (mismatch) CIGAR operators into the ambiguous "
        "'M' operator, for compatibility with tools that do not understand extended CIGAR. By "
        "default, extended CIGAR ('=' / 'X') is used, which lets NM be computed from the CIGAR "
        "alone."
    );
    basic_cigar.option->group( group );
}

// =================================================================================================
//     Settings factory
// =================================================================================================

spear::align::Wfa2Settings Wfa2Options::make_settings() const
{
    spear::align::Wfa2Settings cfg;
    cfg.scope      = spear::align::Wfa2Settings::Scope::Cigar;
    cfg.mismatch   = mismatch_penalty.value;
    cfg.gap_open   = gap_open_penalty.value;
    cfg.gap_extend = gap_extend_penalty.value;

    // --damage-reach and the two individual options are mutually exclusive (enforced by CLI11
    // via excludes() in add_wfa2_opts_to_app()), so at most one of these branches is ever the
    // one the user actually set; the other stays at its unset default of 0.
    if( damage_reach.was_set() ) {
        cfg.damage_ct_5p_reach = damage_reach.value;
        cfg.damage_ga_3p_reach = damage_reach.value;
    } else {
        cfg.damage_ct_5p_reach = damage_ct_5p_reach.value;
        cfg.damage_ga_3p_reach = damage_ga_3p_reach.value;
    }

    cfg.max_steps          = max_steps.value;
    cfg.use_extended_cigar = !basic_cigar.value;
    return cfg;
}
