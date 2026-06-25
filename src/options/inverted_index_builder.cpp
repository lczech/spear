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

#include "options/inverted_index_builder.hpp"

#include "genesis/util/core/exception.hpp"

#include "tools/misc.hpp"

#include <vector>

// =================================================================================================
//      Setup Functions
// =================================================================================================

void InvertedIndexBuilderOptions::add_build_opts_to_app( CLI::App* sub, std::string const& group )
{
    // Correct setup check.
    internal_check(
        !pending_capacity.option && !max_postings_per_term.option && !position_bits.option,
        "Cannot use the same InvertedIndexBuilderOptions object multiple times."
    );

    // Pending capacity.
    pending_capacity = sub->add_option(
        "--pending-capacity",
        pending_capacity.value,
        "Number of positions buffered per k-mer before they are flushed into compressed index "
        "storage. Higher values reduce flush overhead at the cost of more memory."
    );
    pending_capacity.option->check( CLI::NonNegativeNumber );
    pending_capacity.option->group( group );

    // Max occurrences per k-mer (hard cap on stored postings).
    max_postings_per_term = sub->add_option(
        "--max-occurrences-per-kmer",
        max_postings_per_term.value,
        "Maximum number of genome positions to keep per k-mer. If a k-mer occurs more often "
        "than this threshold, all its occurrences are discarded from the index. "
        "The default value of 0 means no limit."
    );
    max_postings_per_term.option->group( group );

    // Position bit width. Hidden, as 32-bit suffices for the vast majority of use cases.
    position_bits = sub->add_option(
        "--position-bits",
        position_bits.value,
        "Number of bits used to store genome bin positions in the index, either 32 or 64. "
        "32 bits suffice for reference genomes of up to about 550 Gbp (with the default "
        "genome bin width); use 64 only for substantially larger references."
    );
    position_bits.option->check( CLI::IsMember( std::vector<std::size_t>{ 32, 64 } ));
    position_bits.option->group( "" );
}
