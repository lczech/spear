#ifndef SPEAR_COMMANDS_MAP_INDEX_H_
#define SPEAR_COMMANDS_MAP_INDEX_H_

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

#include "options/file_output.hpp"
#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

struct MapIndexOptions
{
    // Input
    CliOption<std::string> fasta = "";

    // Settings
    CliOption<size_t>      k = 16;
    CliOption<bool>        canonical = false;
    CliOption<size_t>      max_occurrences_per_kmer = 0;
    CliOption<size_t>      genome_bin_width = 128;

    // Performance
    CliOption<size_t>      pending_capacity = 16;
    CliOption<size_t>      position_bits = 32;

    // Output
    FileOutputOptions      output;
};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_map_index( CLI::App& app );
void run_map_index( MapIndexOptions const& options );

#endif // include guard
