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

#include "commands/tools/version.hpp"
#include "options/global.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/util/core/options.hpp"
#include "genesis/util/text/base64.hpp"

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

// =================================================================================================
//      Setup
// =================================================================================================

// Local forward declaraction.
void run_ee();

void setup_version( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<VersionOptions>();
    auto sub_version = app.add_subcommand(
        "version",
        "Extended version information about spear."
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub_version->callback( [ options ]() {
        run_version( *options );
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_version( VersionOptions const& options )
{
    (void) options;

    LOG_BOLD << spear_header();
    LOG_BOLD;
    LOG_BOLD << "spear version: " << spear_version();
    LOG_BOLD;
    LOG_BOLD << genesis::util::core::info_print_compiler();
    LOG_BOLD << genesis::util::core::info_print_hardware();
    LOG_BOLD;
    LOG_BOLD << "For citation information, call  `spear tools citation`";
    LOG_BOLD << "For license information, call  `spear tools license`";
    LOG_BOLD;
    LOG_BOLD << spear_title();
}
