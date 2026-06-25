#ifndef SPEAR_TOOLS_CLI_OPTION_H_
#define SPEAR_TOOLS_CLI_OPTION_H_

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

#include <string>
#include <type_traits>

// =================================================================================================
//      CLI11 Option Helper
// =================================================================================================

/**
 * @brief Helper that encapsulates an option for the command line interface,
 * storing its value and the CLI11 object used to set that value from the command line.
 *
 * The typical usage pattern is:
 *
 *  1. Declare a member: `CliOption<std::string> out_dir = ".";`
 *  2. Register with CLI11: `out_dir = sub->add_option( "--out-dir", out_dir.value, "..." );`
 *  3. Read the value after parsing: use `out_dir` directly (implicit conversion), or `out_dir.value`.
 *  4. Check if explicitly provided by the user of the program: `out_dir.was_set()`.
 *  5. Chain CLI11 calls via the pointer: `out_dir.option->needs( other_opt );`
 */
template<typename T>
struct CliOption
{
    CliOption() = default;

    /**
     * @brief Construct with an initial value, used as the default before CLI parsing.
     */
    CliOption( T const& val )
        : value( val )
    {}

    /**
     * @brief Construct from a string literal; only available when `T == std::string`.
     */
    template<typename U = T, typename = std::enable_if_t<std::is_same_v<U, std::string>>>
    CliOption( char const* val )
        : value( val )
    {}

    /**
     * @brief Store the CLI11 option pointer returned by `sub->add_option(...)` or `add_flag(...)`.
     *
     * Enables the compact registration idiom: `opt = sub->add_option( name, opt.value, desc );`
     */
    CliOption& operator=( CLI::Option* opt )
    {
        option = opt;
        return *this;
    }

    /**
     * @brief Assign the stored value directly, bypassing the CLI11 option.
     *
     * @note Avoid using the integer literal `0` for `CliOption<bool>` or numeric types: `0` is
     * a valid null pointer constant in C++, so it may resolve to the `CLI::Option*` overload
     * above instead. Use typed literals (`false`, `0u`, etc.) to be unambiguous.
     */
    CliOption& operator=( T const& val )
    {
        value = val;
        return *this;
    }

    /**
     * @brief Assign from a string literal; only available when `T == std::string`.
     */
    template<typename U = T, typename = std::enable_if_t<std::is_same_v<U, std::string>>>
    CliOption& operator=( char const* val )
    {
        value = val;
        return *this;
    }

    /**
     * @brief Implicit conversion to `T const&` for transparent use wherever the value type is expected.
     */
    operator T const&() const
    {
        return value;
    }

    /**
     * @brief Returns true if the option was explicitly provided on the command line,
     * as opposed to just using the default value.
     */
    bool was_set() const
    {
        return option && option->count() > 0;
    }

    /** @brief The option's current value, read by CLI11 and by the rest of the program. */
    T            value  = {};

    /** @brief The CLI11 option object; null until the option is registered via `add_option/add_flag`. */
    CLI::Option* option = nullptr;
};

#endif // include guard
