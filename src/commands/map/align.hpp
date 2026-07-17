#ifndef SPEAR_COMMANDS_MAP_ALIGN_H_
#define SPEAR_COMMANDS_MAP_ALIGN_H_

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

#include "options/bam_writer_options.hpp"
#include "options/file_input.hpp"
#include "options/file_output.hpp"
#include "options/inverted_index_reader.hpp"
#include "options/kmer_seeding_options.hpp"
#include "options/wfa2_options.hpp"

#include "tools/cli_option.hpp"

#include <string>
#include <vector>

// =================================================================================================
//      Options
// =================================================================================================

struct MapAlignOptions
{
    // Input
    InvertedIndexReaderOptions index_reader;
    FileInputOptions           reads;

    // Settings
    KmerSeedingOptions kmer_seeding;
    Wfa2Options        wfa2;

    // Output
    FileOutputOptions output;
    BamWriterOptions   bam_writer;

    // Additionally write the diagnostic TSV (seeding-stage information that has no home in
    // BAM, e.g. num_seeds_found/kept, strand_tie) alongside the BAM output. Off by default;
    // BAM is the primary output.
    CliOption<bool> emit_tsv = false;
};

// =================================================================================================
//      Functions
// =================================================================================================

void setup_map_align( CLI::App& app );
void run_map_align( MapAlignOptions const& options );

#endif // include guard
