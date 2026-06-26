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

#include "commands/map/format.hpp"

#include "spear/reference/functions.hpp"

#include "genesis/util/core/fs.hpp"
#include "genesis/util/format/json/document.hpp"

#include <cstdint>
#include <string>

// =================================================================================================
//     MapIndexMetadata JSON I/O
// =================================================================================================

genesis::util::format::JsonDocument metadata_to_json( MapIndexMetadata const& meta )
{
    using namespace genesis::util::format;

    auto doc = JsonDocument::object();
    doc["index_file"]               = JsonDocument::string( meta.sidx_path );
    doc["k"]                        = JsonDocument::number_unsigned( meta.k );
    doc["canonical"]                = JsonDocument::boolean( meta.canonical );
    doc["genome_bin_width"]         = JsonDocument::number_unsigned( meta.genome_bin_width );
    doc["max_occurrences_per_kmer"] = JsonDocument::number_unsigned( meta.max_occurrences_per_kmer );
    doc["position_bits"]            = JsonDocument::number_unsigned( meta.position_bits );
    doc["total_genome_length"]      = JsonDocument::number_unsigned( meta.total_genome_length );
    return doc;
}

MapIndexMetadata metadata_from_json( genesis::util::format::JsonDocument const& doc )
{
    auto const& obj = doc.get_object();
    MapIndexMetadata meta;
    meta.sidx_path                = obj.at( "index_file" ).get_string();
    meta.k                        = obj.at( "k" ).get_number_unsigned_t<std::uint8_t>();
    meta.canonical                = obj.at( "canonical" ).get_boolean();
    meta.genome_bin_width         = obj.at( "genome_bin_width" ).get_number_unsigned();
    meta.max_occurrences_per_kmer = obj.at( "max_occurrences_per_kmer" ).get_number_unsigned();
    meta.position_bits            = obj.at( "position_bits" ).get_number_unsigned();
    meta.total_genome_length      = obj.at( "total_genome_length" ).get_number_unsigned();
    return meta;
}

// =================================================================================================
//     MapIndexManifest JSON I/O
// =================================================================================================

genesis::util::format::JsonDocument manifest_to_json( MapIndexManifest const& manifest )
{
    using namespace genesis::util::format;
    using namespace spear::reference;

    auto doc         = JsonDocument::object();
    doc["index"]     = metadata_to_json( manifest.metadata );
    doc["reference"] = collection_to_json( manifest.reference );
    return doc;
}

MapIndexManifest manifest_from_json( genesis::util::format::JsonDocument const& doc )
{
    using namespace spear::reference;

    MapIndexManifest manifest;
    auto const& obj    = doc.get_object();
    manifest.metadata  = metadata_from_json( obj.at( "index" ));
    manifest.reference = collection_from_json( obj.at( "reference" ));
    return manifest;
}

// =================================================================================================
//     Path Resolution
// =================================================================================================

std::string resolve_sidx_path( std::string const& json_path, std::string const& sidx_filename )
{
    // Absolute path: use as-is.
    if( !sidx_filename.empty() && sidx_filename[0] == '/' ) {
        return sidx_filename;
    }
    auto const dir = genesis::util::core::file_path( json_path );
    return dir.empty() ? sidx_filename : dir + "/" + sidx_filename;
}
