#ifndef SPEAR_OPTIONS_KMER_SEEDING_H_
#define SPEAR_OPTIONS_KMER_SEEDING_H_

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

#include "spear/seeding/kmer_seeding.hpp"

#include "tools/cli_option.hpp"

#include <cstddef>
#include <string>

// =================================================================================================
//     KmerSeedingOptions
// =================================================================================================

/**
 * @brief CLI binding for KmerSeeding configuration.
 *
 * Pure CLI wrapper: registers the k-mer seeding options with a CLI11 subcommand and converts
 * them to a KmerSeeding::Config. Use make_engine() to construct the runtime engine.
 *
 * Embed one instance in a command's options struct and call add_kmer_seeding_opts_to_app() once.
 */
class KmerSeedingOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    KmerSeedingOptions()  = default;
    ~KmerSeedingOptions() = default;

    KmerSeedingOptions( KmerSeedingOptions const& ) = default;
    KmerSeedingOptions( KmerSeedingOptions&& )      = default;

    KmerSeedingOptions& operator=( KmerSeedingOptions const& ) = default;
    KmerSeedingOptions& operator=( KmerSeedingOptions&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup
    // -------------------------------------------------------------------------

    /**
     * @brief Register all k-mer seeding options with @p sub.
     *
     * Must not be called more than once on the same instance.
     */
    void add_kmer_seeding_opts_to_app( CLI::App* sub, std::string const& group = "K-mer Seeding" );

    // -------------------------------------------------------------------------
    //     Engine factory
    // -------------------------------------------------------------------------

    /**
     * @brief Construct a KmerSeeding engine from the current option values.
     *
     * @param index      Opened InvertedIndex; must outlive the returned engine.
     * @param bin_width  Genome bin width used when the index was built.
     */
    template<typename PositionT>
    spear::seeding::KmerSeeding<PositionT> make_engine(
        spear::inverted_index::InvertedIndex<PositionT> const& index,
        std::size_t bin_width
    ) const {
        typename spear::seeding::KmerSeeding<PositionT>::Config cfg;
        cfg.min_hit_count            = min_hit_count.value;
        cfg.min_hit_fraction         = min_hit_fraction.value;
        cfg.max_occurrences_per_kmer = max_occurrences_per_kmer.value;
        return spear::seeding::KmerSeeding<PositionT>( cfg, index, bin_width );
    }

    // -------------------------------------------------------------------------
    //     Option members
    // -------------------------------------------------------------------------

    // Minimum number of distinct k-mer hits within a window; 0 = auto-derive per read.
    CliOption<std::size_t> min_hit_count          = 0;

    // Fraction of a read's k-mers needed when auto-deriving min_hit_count.
    CliOption<double>      min_hit_fraction       = 0.0;

    // Soft cap on per-k-mer occurrences; 0 = no cap.
    CliOption<std::size_t> max_occurrences_per_kmer = 0;

};

#endif // SPEAR_OPTIONS_KMER_SEEDING_H_
