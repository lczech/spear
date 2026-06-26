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

#include "options/inverted_index_reader.hpp"

#include "tools/misc.hpp"

#include <vector>

// =================================================================================================
//     Setup
// =================================================================================================

void InvertedIndexReaderOptions::add_reader_opts_to_app(
    CLI::App* sub, std::string const& group
) {
    internal_check(
        !index_file.option && !load_mode.option,
        "Cannot use the same InvertedIndexReaderOptions object multiple times."
    );

    // Index file path.
    index_file = sub->add_option(
        "--index-file",
        index_file.value,
        "Path to the JSON manifest of the index, as produced by `spear map index`. "
        "The binary index data (.sidx) is located automatically from the path recorded "
        "in the manifest."
    );
    index_file.option->check( CLI::ExistingFile );
    index_file.option->required();
    index_file.option->group( group );

    // Load mode. Hidden: most users never need to change this.
    load_mode = sub->add_option(
        "--index-load-mode",
        load_mode.value,
        "How to access the index's binary posting-list data: 'load-all' reads the entire blob "
        "into memory at startup (recommended); 'pread' keeps it on disk and issues one pread() "
        "syscall per k-mer lookup, suitable when the blob does not fit in memory, "
        "at the exense of performance."
    );
    load_mode.option->check( CLI::IsMember( std::vector<std::string>{ "load-all", "pread" } ));
    load_mode.option->group( "" );
}
