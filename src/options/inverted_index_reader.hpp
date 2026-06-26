#ifndef SPEAR_OPTIONS_INVERTED_INDEX_READER_H_
#define SPEAR_OPTIONS_INVERTED_INDEX_READER_H_

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

#include <string>

// =================================================================================================
//     InvertedIndexReaderOptions
// =================================================================================================

/**
 * @brief Options helper for commands that open an existing inverted index.
 *
 * Wraps the `--index-file` path option and the hidden `--index-load-mode` choice.
 * Embed one instance in a command's options struct and call add_reader_opts_to_app() once.
 */
class InvertedIndexReaderOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    InvertedIndexReaderOptions()  = default;
    ~InvertedIndexReaderOptions() = default;

    InvertedIndexReaderOptions( InvertedIndexReaderOptions const& ) = default;
    InvertedIndexReaderOptions( InvertedIndexReaderOptions&& )      = default;

    InvertedIndexReaderOptions& operator=( InvertedIndexReaderOptions const& ) = default;
    InvertedIndexReaderOptions& operator=( InvertedIndexReaderOptions&& )      = default;

    // -------------------------------------------------------------------------
    //     Setup
    // -------------------------------------------------------------------------

    /**
     * @brief Register `--index-file` and the hidden `--index-load-mode` with @p sub.
     *
     * Must not be called more than once on the same instance.
     */
    void add_reader_opts_to_app( CLI::App* sub, std::string const& group = "Input" );

    // -------------------------------------------------------------------------
    //     Open Mode Helper
    // -------------------------------------------------------------------------

    /**
     * @brief Return true if the parsed load mode is `pread` (keep blob on disk).
     *
     * Use this to select InvertedIndex<PositionT>::OpenMode at call sites where the template
     * parameter is already known:
     * @code
     *   auto mode = opts.is_pread()
     *       ? InvertedIndex<PositionT>::OpenMode::kPread
     *       : InvertedIndex<PositionT>::OpenMode::kLoadAll;
     * @endcode
     */
    bool is_pread() const
    {
        return load_mode.value == "pread";
    }

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

    // Path to the JSON manifest file produced by `spear map index`.
    CliOption<std::string> index_file = "";

    // Load mode for the binary posting-list blob: "load-all" (default) or "pread".
    CliOption<std::string> load_mode  = "load-all";
};

#endif // SPEAR_OPTIONS_INVERTED_INDEX_READER_H_
