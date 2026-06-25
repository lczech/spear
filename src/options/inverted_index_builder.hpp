#ifndef SPEAR_OPTIONS_INVERTED_INDEX_BUILDER_H_
#define SPEAR_OPTIONS_INVERTED_INDEX_BUILDER_H_

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

#include "tools/cli_option.hpp"

#include <cstddef>
#include <string>

// =================================================================================================
//      Inverted Index Builder Options
// =================================================================================================

/**
 * @brief Options helper class for building an inverted index via InvertedIndexBuilder.
 *
 * Wraps the user-controllable fields of InvertedIndexBuilder::Config and the position-bit-width
 * choice that selects the PositionT template parameter. One instance is meant to be embedded in
 * a command's options struct and registered with a single add_build_opts_to_app() call.
 *
 * After CLI parsing, call make_builder_config<Builder>() to obtain the fully-typed
 * InvertedIndexBuilder::Config, and branch on position_bits to select the Builder instantiation.
 */
class InvertedIndexBuilderOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    InvertedIndexBuilderOptions()  = default;
    ~InvertedIndexBuilderOptions() = default;

    InvertedIndexBuilderOptions( InvertedIndexBuilderOptions const& ) = default;
    InvertedIndexBuilderOptions( InvertedIndexBuilderOptions&& )      = default;

    InvertedIndexBuilderOptions& operator=( InvertedIndexBuilderOptions const& ) = default;
    InvertedIndexBuilderOptions& operator=( InvertedIndexBuilderOptions&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Register all builder options with the given CLI11 sub-app.
     *
     * Adds `--pending-capacity` and `--max-occurrences-per-kmer` (both visible, in @p group) and
     * `--position-bits` (hidden). Must not be called more than once on the same instance.
     */
    void add_build_opts_to_app( CLI::App* sub, std::string const& group = "Index Settings" );

    // -------------------------------------------------------------------------
    //     Builder Config Helper
    // -------------------------------------------------------------------------

    /**
     * @brief Construct an InvertedIndexBuilder::Config from the parsed option values.
     *
     * @tparam Builder  The fully-specialised InvertedIndexBuilder type to use.
     *                  position_bits is not reflected here; the caller is responsible for
     *                  branching on position_bits.value to select the correct Builder type.
     */
    template<typename Builder>
    typename Builder::Config make_builder_config() const
    {
        typename Builder::Config cfg;
        cfg.pending_capacity     = pending_capacity.value;
        cfg.max_postings_per_term = static_cast<typename Builder::stored_count_type>(
            max_postings_per_term.value
        );
        return cfg;
    }

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

public:

    // Number of positions buffered per k-mer before flushing into compressed storage.
    CliOption<std::size_t> pending_capacity    = std::size_t{ 16 };

    // Hard cap on unique positions per k-mer; 0 means no cap.
    CliOption<std::size_t> max_postings_per_term = std::size_t{ 0 };

    // Bit width of the position type stored in the index: 32 or 64.
    CliOption<std::size_t> position_bits       = std::size_t{ 32 };

};

#endif // include guard
