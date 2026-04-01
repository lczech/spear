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

#include "commands/tools/version.hpp"
#include "options/global.hpp"
#include "tools/references.hpp"
#include "tools/version.hpp"

#include "genesis/util/core/options.hpp"
#include "genesis/util/text/base64.hpp"

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

// =================================================================================================
//      Setup
// =================================================================================================

// Local forward declaraction.
void run_ee();

void setup_version( CLI::App& app )
{
    // Create the options and subcommand objects.
    auto options = std::make_shared<VersionOptions>();
    auto sub_version = app.add_subcommand(
        "version",
        "Extended version information about spear."
    );

    // Set the run function as callback to be called when this subcommand is issued.
    // Hand over the options by copy, so that their shared ptr stays alive in the lambda.
    sub_version->callback( [ options ]() {
        run_version( *options );
    });

    // Same for the ee command, but this time without options at all.
    auto sub_ee = app.add_subcommand(
        genesis::util::text::base64_decode_string(
            "cXVvdGU="
        ),
        genesis::util::text::base64_decode_string(
            "R2V0IHRoZSB3aXNkb20gb2YgdGhlIGFuY2llbnQgR3JlZWsgcGhpbG9zb3BoZXJzLg=="
        )
    )->group( "" );
    sub_ee->callback( []() {
        run_ee();
    });
}

// =================================================================================================
//      Run
// =================================================================================================

void run_version( VersionOptions const& options )
{
    (void) options;

    LOG_BOLD << spear_header();
    LOG_BOLD;
    LOG_BOLD << "spear version: " << spear_version();
    LOG_BOLD;
    LOG_BOLD << genesis::util::core::info_print_compiler();
    LOG_BOLD << genesis::util::core::info_print_hardware();
    LOG_BOLD;
    LOG_BOLD << "For citation information, call  `spear tools citation`";
    LOG_BOLD << "For license information, call  `spear tools license`";
    LOG_BOLD;
    LOG_BOLD << spear_title();
}

void run_ee()
{
    std::vector<std::string> qs = {
        "Ikkga25vdyB0aGF0IEkga25vdyBub3RoaW5nLiIg4oCTIFNvY3JhdGVz", "IkFuIHVuZXhhbWluZWQgbGlmZSBpcy"
        "Bub3Qgd29ydGggbGl2aW5nLiIg4oCTIFNvY3JhdGVz", "IlRoZSBvbmx5IHRydWUgd2lzZG9tIGlzIGluIGtub3dp"
        "bmcgeW91IGtub3cgbm90aGluZy4iIOKAkyBTb2NyYXRlcw==", "Ildpc2UgbWVuIHNwZWFrIGJlY2F1c2UgdGhleS"
        "BoYXZlIHNvbWV0aGluZyB0byBzYXk7IGZvb2xzIGJlY2F1c2UgdGhleSBoYXZlIHRvIHNheSBzb21ldGhpbmcuIiDi"
        "gJMgUGxhdG8=","IldlIGFyZSB3aGF0IHdlIHJlcGVhdGVkbHkgZG8uIEV4Y2VsbGVuY2UsIHRoZW4sIGlzIG5vdCB"
        "hbiBhY3QsIGJ1dCBhIGhhYml0LiIg4oCTIEFyaXN0b3RsZQ==", "IlRoZXJlIGlzIG5vdGhpbmcgcGVybWFuZW50I"
        "GV4Y2VwdCBjaGFuZ2UuIiDigJMgSGVyYWNsaXR1cw==", "Ik5vIG1hbiBldmVyIHN0ZXBzIGluIHRoZSBzYW1lIHJ"
        "pdmVyIHR3aWNlLiIg4oCTIEhlcmFjbGl0dXM=", "IkNoYXJhY3RlciBpcyBkZXN0aW55LiIg4oCTIEhlcmFjbGl0d"
        "XM=", "IlRoZSBtaW5kIGlzIG5vdCBhIHZlc3NlbCB0byBiZSBmaWxsZWQgYnV0IGEgZmlyZSB0byBiZSBraW5kbGV"
        "kLiIg4oCTIFBsdXRhcmNo", "Iktub3cgaG93IHRvIGxpc3RlbiBhbmQgeW91IHdpbGwgcHJvZml0IGV2ZW4gZnJvb"
        "SB0aG9zZSB3aG8gdGFsayBiYWRseS4iIOKAkyBQbHV0YXJjaA==", "IldlYWx0aCBjb25zaXN0cyBub3QgaW4gaGF"
        "2aW5nIGdyZWF0IHBvc3Nlc3Npb25zLCBidXQgaW4gaGF2aW5nIGZldyB3YW50cy4iIOKAkyBFcGljdGV0dXM=", "I"
        "k5vIG1hbiBpcyBmcmVlIHdobyBpcyBub3QgYSBtYXN0ZXIgb2YgaGltc2VsZi4iIOKAkyBFcGljdGV0dXM="
    };

    auto const seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine engine( seed );
    std::shuffle( qs.begin(), qs.end(), engine );
    LOG_BOLD << genesis::util::text::base64_decode_string( qs[0] );
}
