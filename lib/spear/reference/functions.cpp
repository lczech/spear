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

#include "spear/reference/functions.hpp"

#include "genesis/sequence/format/fasta_reader.hpp"
#include "genesis/sequence/format/fastx_input_stream.hpp"
#include "genesis/sequence/function/dict.hpp"
#include "genesis/util/core/logging.hpp"
#include "genesis/util/io/input_source.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace spear {
namespace reference {

// =================================================================================================
//     Builder
// =================================================================================================

ReferenceCollection build_collection(
    std::vector<std::string> const& fasta_paths,
    bool store_sites
) {
    using namespace genesis::sequence;
    using namespace genesis::util;

    ReferenceCollection result;

    for( auto const& fasta_path : fasta_paths ) {
        if( store_sites ) {
            // Always read full FASTA when sequence data is requested.
            for( auto& seq : FastaInputStream( io::from_file( fasta_path ))) {
                result.add( fasta_path, std::move( seq ) );
            }
            continue;
        }

        // Try FAI first: fastest path, no sequence data needed at all.
        auto const fai_path = genesis::sequence::find_sequence_fai( fasta_path );
        if( !fai_path.empty() ) {
            auto const dict = read_sequence_fai( io::from_file( fai_path ) );
            for( auto const& entry : dict ) {
                result.add( fasta_path, entry.name, static_cast<std::uint64_t>( entry.length ) );
            }
            continue;
        }

        // Try sequence dictionary next.
        auto const dict_path = genesis::sequence::find_sequence_dict( fasta_path );
        if( !dict_path.empty() ) {
            auto const dict = read_sequence_dict( io::from_file( dict_path ) );
            for( auto const& entry : dict ) {
                result.add( fasta_path, entry.name, static_cast<std::uint64_t>( entry.length ) );
            }
            continue;
        }

        // Fall back to a full FASTA scan; warn so the user can create a FAI for next time.
        LOG_WARN
            << "No FAI index or sequence dictionary found for '" << fasta_path << "'. "
            << "Scanning the full FASTA file to obtain sequence names and lengths. "
            << "Run `samtools faidx " << fasta_path << "` to create an index "
            << "and avoid this scan on subsequent runs."
        ;
        for( auto& seq : FastaInputStream( io::from_file( fasta_path ))) {
            result.add( fasta_path, seq.label(), static_cast<std::uint64_t>( seq.length() ) );
        }
    }

    return result;
}

// =================================================================================================
//     JSON I/O
// =================================================================================================

genesis::util::format::JsonDocument to_json( ReferenceCollection const& collection )
{
    using namespace genesis::util::format;

    // Fasta files list.
    auto files_arr = JsonDocument::array();
    for( auto const& f : collection.fasta_files() ) {
        files_arr.push_back( JsonDocument::string( f ) );
    }

    // Sequences list. Sites are not stored in JSON (too large); file_index and length are enough
    // to reconstruct the collection, with global offsets re-derived from the prefix sums.
    auto seqs_arr = JsonDocument::array();
    for( std::size_t i = 0; i < collection.size(); ++i ) {
        auto const& s = collection[i];
        auto entry = JsonDocument::object();
        entry["header"]     = JsonDocument::string( s.header );
        entry["file_index"] = JsonDocument::number_unsigned( s.file_index );
        entry["length"]     = JsonDocument::number_unsigned( s.length );
        seqs_arr.push_back( std::move( entry ) );
    }

    auto doc = JsonDocument::object();
    doc["fasta_files"] = std::move( files_arr );
    doc["sequences"]   = std::move( seqs_arr );
    return doc;
}

ReferenceCollection from_json( genesis::util::format::JsonDocument const& doc )
{
    auto const& obj = doc.get_object();

    // Read the flat fasta_files list so per-sequence file_index can be resolved.
    std::vector<std::string> files;
    for( auto const& f : obj.at( "fasta_files" ).get_array() ) {
        files.push_back( f.get_string() );
    }

    // Reconstruct sequences via the public add() interface; lengths drive the prefix-sum offsets.
    ReferenceCollection result;
    for( auto const& entry : obj.at( "sequences" ).get_array() ) {
        auto const& e          = entry.get_object();
        auto const  header     = e.at( "header" ).get_string();
        auto const  file_index = e.at( "file_index" ).get_number_unsigned();
        auto const  length     = e.at( "length" ).get_number_unsigned();
        result.add( files[file_index], header, length );
    }

    return result;
}

} // namespace reference
} // namespace spear
