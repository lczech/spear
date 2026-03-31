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

// #include "commands/commands.hpp"
#include "commands/tools.hpp"

#include "options/global.hpp"

#include "tools/cli_formatter.hpp"
#include "tools/cli_setup.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/util/core/exception.hpp"
#include "genesis/util/core/logging.hpp"

#include <memory>
#include <sstream>

// =================================================================================================
//      Main Program
// =================================================================================================

int main( int argc, char** argv )
{
    // -------------------------------------------------------------------------
    //     Logging
    // -------------------------------------------------------------------------

    // Activate logging.
    genesis::util::core::Logging::log_to_stdout();
    genesis::util::core::Logging::details.level = false;
    // genesis::util::core::Logging::details.time = true;

    // -------------------------------------------------------------------------
    //     App Setup
    // -------------------------------------------------------------------------

    // Set up the main CLI app.
    CLI::App app{};
    app.name( "spear" );
    app.description( spear_header() );
    app.footer( "\n" + spear_title() );

    // We use a custom formatter for the help messages that looks nicer.
    app.formatter( std::make_shared<SpearFormatter>() );
    app.get_formatter()->enable_description_formatting(false);

    // We don't like short options in spear. Reset the help option.
    // This is inherited by subcommands automatically.
    app.set_help_flag( "--help", "Print this help message and exit." );
    app.set_version_flag( "--version", spear_version(), "Print the spear version and exit." );

    // We want all options to capture their default values, so that we see them in the help.
    // No idea why CLI does not do that all the time. Hopefully, this is the correct way to solve this.
    app.option_defaults()->always_capture_default( true );

    // spear expects exactly one subcommand.
    app.require_subcommand( 0, 1 );

    // -------------------------------------------------------------------------
    //     Subcommand Setup
    // -------------------------------------------------------------------------

    // Init global options. This ensures that all subcommands that need the global options
    // have them initialized to proper values that they can use (e.g., for the help output).
    global_options.initialize( argc, argv );

    // Set up all subcommands.
    // setup_commands( app );
    setup_tools( app );

    // -------------------------------------------------------------------------
    //     Final Checks and Steps
    // -------------------------------------------------------------------------

    // General checks before running.
    // These are mainly to support the development, to avoid bugs and mistakes.
    check_all_citations();
    check_unique_command_names( app );
    check_subcommand_names( app );

    // Final steps before we can run.
    // Make sure that we have all defaults captured,
    // so that they can be used in the header print of each command.
    fix_cli_default_values( app );

    // -------------------------------------------------------------------------
    //     Go Go Go
    // -------------------------------------------------------------------------

    int exit_code = 0;
    try {

        // Do the parsing of the command line arguments.
        // Because we use callback functions for the commands,
        // this is also where the actual commands and functions are executed.
        app.parse( argc, argv );

    } catch( CLI::ParseError const& error ) {

        // Capture the app exit message, so that we can print it to our log.
        std::stringstream ss;
        exit_code = app.exit( error, ss, ss );
        auto const message = ss.str();
        if( ! message.empty() ) {
            LOG_BOLD << message;
        }

    } catch( genesis::util::core::ExistingFileError const& error ) {

        // Special case for existing files: This is very common, and we want a nice and useful
        // error messsage for this one!
        std::string message
            = "Output file '" + error.filename() + "' already exists. If you want to allow "
            + "overwriting of existing output files, use "
            + global_options.opt_allow_file_overwriting.option->get_name()
        ;
        LOG_BOLD << message;
        LOG_BOLD;
        throw genesis::util::core::ExistingFileError( message, error.filename() );

    } catch( std::exception const& error ) {

        // More general error. We also catch this in order to make sure it is logged everywhere,
        // but then re-throw to make the program fail.
        LOG_BOLD << "Error: " << error.what();
        LOG_BOLD;
        throw;

    } catch (...) {

        // Capture everything else. Should not happen, but why not.
        LOG_BOLD << "Error: Unknown error.";
        LOG_BOLD;
        throw;
    }

    // Close all logging, and return the exit code.
    genesis::util::core::Logging::clear();
    return exit_code;
}
