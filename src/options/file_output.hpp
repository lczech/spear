#ifndef SPEAR_OPTIONS_FILE_OUTPUT_H_
#define SPEAR_OPTIONS_FILE_OUTPUT_H_

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

#include "genesis/util/io/output_target.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// =================================================================================================
//      File Output Options
// =================================================================================================

/**
 * @brief Options helper class to set up consistent file output.
 *
 * One instance of this class is meant to be used for a file or set of files that shall go into
 * the same output directory. For complex commands that produce several sets of files, it is
 * recommended to use multiple instances of this class, so that the user can provide separate
 * output directories for each set.
 *
 * Furthermore, the compress option also affects all files that are added via this class.
 * Hence, if you want to offer to compress only a particular file, but not all, again, use
 * separate instances of this class.
 */
class FileOutputOptions
{
public:

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    FileOutputOptions()  = default;
    virtual ~FileOutputOptions() = default;

    FileOutputOptions( FileOutputOptions const& other ) = default;
    FileOutputOptions( FileOutputOptions&& )            = default;

    FileOutputOptions& operator= ( FileOutputOptions const& other ) = default;
    FileOutputOptions& operator= ( FileOutputOptions&& )            = default;

    // -------------------------------------------------------------------------
    //     File Type/Name Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Set the name of the options added via this class.
     *
     * This name is used as part of the command line option name, for example `--optionname-out-dir`.
     * It has to be set prior to calling any of the `add...` functions of this class, so that
     * those can actually use the optionname.
     *
     * Default is empty, creating options of the form `--out-dir`, without any extra option name.
     */
    void set_optionname( std::string const& optionname );

    /**
     * @brief Set the group for the options added via this class.
     *
     * It has to be set prior to calling any of the `add...` functions of this class, so that
     * those can actually use the group name.
     *
     * Default is "Output".
     */
    void set_group( std::string const& group );

    // -------------------------------------------------------------------------
    //     Extra Setup Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Shortcut to set three options at the same time, out-dir, prefix, and suffix.
     */
    void add_default_output_opts_to_app(
        CLI::App* sub,
        std::string const& output_dir_initial_value = ".",
        std::string const& prefix_initial_value = "",
        std::string const& suffix_initial_value = ""
    );

    /**
     * @brief Add a named output dir "--name-out-dir", or "--out-dir" if name is empty.
     */
    CLI::Option* add_output_dir_opt_to_app(
        CLI::App* sub,
        std::string const& initial_value = "."
    );

    /**
     * @brief Add a command to set the file prefix for this file name, using "--name-file-prefix",
     * or "--file-prefix" if name is empty.
     */
    CLI::Option* add_file_prefix_opt_to_app(
        CLI::App* sub,
        std::string const& initial_value = ""
    );

    /**
     * @brief Add a command to set the file suffix for this file name, using "--name-file-suffix",
     * or "--file-suffix" if name is empty.
     */
    CLI::Option* add_file_suffix_opt_to_app(
        CLI::App* sub,
        std::string const& initial_value = ""
    );

private:

    /**
     * @brief Add either the prefix or the suffix option.
     */
    CLI::Option* add_filefix_opt_(
        CLI::App* sub,
        std::string const& initial_value,
        std::string const& fixname,
        CliOption<std::string>& target
    );

public:

    /**
     * @brief Add an option to produce gzip compressed output, using "--name-compress" if optionname
     * was given in the set_optionname() call, or "--compress" if optionname was empty.
     */
    CLI::Option* add_file_compress_opt_to_app(
        CLI::App* sub
    );

    // -------------------------------------------------------------------------
    //     Accessors
    // -------------------------------------------------------------------------

    std::string const& optionname() const
    {
        return optionname_;
    }

    std::string const& group() const
    {
        return group_;
    }

    // -------------------------------------------------------------------------
    //     Output Filenames
    // -------------------------------------------------------------------------

    std::string get_output_filename(
        std::string const& infix, std::string const& extension, bool with_dir = true
    ) const;

    void check_output_files_nonexistence(
        std::string const& infix, std::string const& extension
    ) const;

    void check_output_files_nonexistence(
        std::vector<std::string> const& infixes, std::string const& extension
    ) const;

    void check_output_files_nonexistence(
        std::string const& infix, std::vector<std::string> const& extensions
    ) const;

    void check_output_files_nonexistence(
        std::vector<std::pair<std::string, std::string>> const& infixes_and_extensions
    ) const;

    // -------------------------------------------------------------------------
    //     Output Targets
    // -------------------------------------------------------------------------

    std::shared_ptr<genesis::util::io::BaseOutputTarget> get_output_target(
        std::string const& infix, std::string const& extension
    ) const;

    // -------------------------------------------------------------------------
    //     Option Members
    // -------------------------------------------------------------------------

public:

    CliOption<std::string> out_dir  = ".";
    CliOption<std::string> prefix;
    CliOption<std::string> suffix;
    CliOption<bool>        compress = false;

private:

    // Internal config; must be set before any add_*_opt_to_app() calls.
    std::string optionname_ = "";
    std::string group_ = "Output";

};

#endif // include guard
