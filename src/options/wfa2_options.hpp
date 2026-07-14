#ifndef SPEAR_OPTIONS_WFA2_H_
#define SPEAR_OPTIONS_WFA2_H_

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

#include "CLI/CLI.hpp"

#include "spear/align/wfa2.hpp"

#include "tools/cli_option.hpp"

#include <string>

// =================================================================================================
//     Wfa2Options
// =================================================================================================

/**
 * @brief CLI binding for Wfa2Aligner configuration.
 *
 * Pure CLI wrapper: registers the WFA2 alignment options with a CLI11 subcommand and converts
 * them to a Wfa2Settings. Unlike KmerSeedingOptions::make_engine(), this stops one level short
 * of a ready-to-use engine: Wfa2Aligner is not thread-safe (settings are fixed at construction,
 * one instance per thread), so callers of make_settings() construct their own Wfa2Aligner at
 * the call site -- typically as a thread_local -- from the returned settings.
 *
 * Embed one instance in a command's options struct and call add_wfa2_opts_to_app() once.
 */
class Wfa2Options
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    Wfa2Options()  = default;
    ~Wfa2Options() = default;

    Wfa2Options( Wfa2Options const& ) = default;
    Wfa2Options( Wfa2Options&& )      = default;

    Wfa2Options& operator=( Wfa2Options const& ) = default;
    Wfa2Options& operator=( Wfa2Options&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup
    // -------------------------------------------------------------------------

    /**
     * @brief Register all WFA2 alignment options with @p sub.
     *
     * Must not be called more than once on the same instance.
     */
    void add_wfa2_opts_to_app( CLI::App* sub, std::string const& group = "Alignment" );

    // -------------------------------------------------------------------------
    //     Settings factory
    // -------------------------------------------------------------------------

    /**
     * @brief Build a Wfa2Settings from the current option values.
     *
     * Always requests Wfa2Settings::Scope::Cigar (not exposed as an option): this is the only
     * scope with a well-defined alignment position (ref_begin requires the backtrace that
     * Scope::Score skips), so it is what every caller needs for real alignment output.
     */
    spear::align::Wfa2Settings make_settings() const;

    // -------------------------------------------------------------------------
    //     Option members
    // -------------------------------------------------------------------------

    // Mismatch penalty. Must be > 0: WFA2 itself requires this (see wfa2.cpp).
    CliOption<int> mismatch_penalty   = 4;

    // Gap-open penalty. May be 0 (unlike the mismatch/gap-extend penalties).
    CliOption<int> gap_open_penalty   = 4;

    // Gap-extension penalty. Must be > 0: WFA2 itself requires this (see wfa2.cpp).
    CliOption<int> gap_extend_penalty = 2;

    // Bases from the 5' end where C→T is treated as a free match (aDNA damage tolerance).
    // 0 disables 5' damage tolerance (default). Mutually exclusive with damage_reach.
    CliOption<int> damage_ct_5p_reach = 0;

    // Bases from the 3' end where G→A is treated as a free match (aDNA damage tolerance).
    // 0 disables 3' damage tolerance (default). Mutually exclusive with damage_reach.
    CliOption<int> damage_ga_3p_reach = 0;

    // Convenience option that sets both damage_ct_5p_reach and damage_ga_3p_reach to the same
    // value; mutually exclusive with setting them individually. 0 (default, unset) means: use
    // damage_ct_5p_reach/damage_ga_3p_reach instead.
    CliOption<int> damage_reach = 0;

    // Safety cap on alignment score (WFA steps) before aborting a pathological alignment.
    CliOption<int> max_steps = 4096;

    // If set, collapse the X/= CIGAR operators into the ambiguous M operator, for compatibility
    // with tools that do not understand extended CIGAR. Off by default (extended CIGAR is used).
    CliOption<bool> basic_cigar = false;

};

#endif // SPEAR_OPTIONS_WFA2_H_
