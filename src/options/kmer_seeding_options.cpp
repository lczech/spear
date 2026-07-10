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

#include "options/kmer_seeding_options.hpp"

#include "tools/misc.hpp"

// =================================================================================================
//     Setup
// =================================================================================================

void KmerSeedingOptions::add_kmer_seeding_opts_to_app(
    CLI::App* sub, std::string const& group
) {
    internal_check(
        !min_hit_count.option && !min_hit_fraction.option && !max_occurrences_per_kmer.option,
        "Cannot use the same KmerSeedingOptions object multiple times."
    );

    // Min hit count.
    min_hit_count = sub->add_option(
        "--min-hit-count",
        min_hit_count.value,
        "Minimum number of distinct k-mer hits required within a window for it to count as a "
        "candidate interval. The default value of 0 derives this automatically per read: "
        "if --min-hit-fraction is also 0 (the default), any window with at least 1 hit is "
        "reported; otherwise round(num_kmers * --min-hit-fraction) is used."
    );
    min_hit_count.option->check( CLI::NonNegativeNumber );
    min_hit_count.option->group( group );

    // Min hit fraction.
    min_hit_fraction = sub->add_option(
        "--min-hit-fraction",
        min_hit_fraction.value,
        "Fraction of a read's k-mers that need to be found within a window for it to count as a "
        "candidate interval. Used to automatically derive --min-hit-count when that option is "
        "left at its default of 0. A value of 0.0 (the default) means use at least 1 hit."
    );
    min_hit_fraction.option->check( CLI::Range( 0.0, 1.0 ));
    min_hit_fraction.option->group( group );

    // Max occurrences per k-mer.
    max_occurrences_per_kmer = sub->add_option(
        "--max-occurrences-per-kmer",
        max_occurrences_per_kmer.value,
        "Soft cap on the number of occurrences of a k-mer: k-mers that occur more often than "
        "this are skipped without decoding, as if their data were absent from the index. Unlike "
        "`spear map index`'s `--max-occurrences-per-kmer` (a hard cap baked into the index at "
        "build time), this can be changed freely per run to experiment with different limits "
        "without rebuilding the index. The default value of 0 means no soft cap is applied."
    );
    max_occurrences_per_kmer.option->check( CLI::NonNegativeNumber );
    max_occurrences_per_kmer.option->group( group );
}
