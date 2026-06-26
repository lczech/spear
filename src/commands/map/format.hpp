#ifndef SPEAR_COMMANDS_MAP_FORMAT_H_
#define SPEAR_COMMANDS_MAP_FORMAT_H_

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

/**
 * @brief Data structures and JSON I/O for the spear map-index manifest file.
 *
 * @file
 * @ingroup formats
 */

#include "spear/reference/collection.hpp"

#include "genesis/util/format/json/document.hpp"

#include <cstdint>
#include <string>

// =================================================================================================
//     MapIndexMetadata
// =================================================================================================

/**
 * @brief Index build parameters for a `spear map index`, stored in the `"index"` sub-object
 * of the JSON manifest file.
 *
 * Written by `spear map index` via metadata_to_json(), read back via metadata_from_json().
 * Lives at the command/format layer rather than inside the inverted-index library, keeping
 * that library free of genomics terminology.
 *
 * See MapIndexManifest for the full manifest structure combining this with a ReferenceCollection.
 */
struct MapIndexMetadata
{
    // Basename of the binary index file (.sidx), relative to the JSON manifest.
    std::string sidx_path;

    // K-mer size used when building the index.
    std::uint8_t k = 16;

    // Whether canonical (strand-collapsed) k-mers were indexed.
    bool canonical = false;

    // Width of the genome bins (in bases) that positions are grouped into.
    std::uint64_t genome_bin_width = 128;

    // Hard cap on unique positions per k-mer baked in at build time; 0 = no cap.
    std::uint64_t max_occurrences_per_kmer = 0;

    // Bit width of stored positions in the binary index: 32 or 64.
    std::uint64_t position_bits = 32;

    // Total length in bases of all reference sequences indexed.
    std::uint64_t total_genome_length = 0;
};

// =================================================================================================
//     MapIndexManifest
// =================================================================================================

/**
 * @brief Full contents of a `spear map index` JSON manifest file.
 *
 * Combines the index build parameters (MapIndexMetadata) with the reference sequence catalog
 * (ReferenceCollection). Written by `spear map index` via manifest_to_json(), read by
 * `spear map locate` (and future commands) via manifest_from_json().
 *
 * The JSON has two top-level sub-objects:
 *   - `"index"`:     MapIndexMetadata, serialised by metadata_to_json().
 *   - `"reference"`: ReferenceCollection, serialised by spear::reference::collection_to_json().
 */
struct MapIndexManifest
{
    // Index build parameters (the "index" sub-object).
    MapIndexMetadata metadata;

    // Reference sequence catalog (the "reference" sub-object).
    spear::reference::ReferenceCollection reference;
};

// =================================================================================================
//     Free Functions
// =================================================================================================

/**
 * @brief Serialize MapIndexMetadata to a JSON object for the `"index"` sub-object.
 */
genesis::util::format::JsonDocument metadata_to_json( MapIndexMetadata const& meta );

/**
 * @brief Deserialize MapIndexMetadata from the `"index"` sub-object of a map-index manifest.
 *
 * @p doc must be the `"index"` sub-object (e.g., `manifest_doc.get_object().at("index")`),
 * not the top-level manifest document.
 */
MapIndexMetadata metadata_from_json( genesis::util::format::JsonDocument const& doc );

/**
 * @brief Serialize a MapIndexManifest to a JSON document with `"index"` and `"reference"`
 * sub-objects, suitable for writing as the map-index manifest file.
 */
genesis::util::format::JsonDocument manifest_to_json( MapIndexManifest const& manifest );

/**
 * @brief Deserialize a MapIndexManifest from the top-level JSON manifest document.
 *
 * @p doc must be the top-level document (with `"index"` and `"reference"` keys),
 * as produced by manifest_to_json().
 */
MapIndexManifest manifest_from_json( genesis::util::format::JsonDocument const& doc );

/**
 * @brief Resolve the `.sidx` binary path relative to the directory of its JSON manifest.
 *
 * If @p sidx_filename is absolute it is returned unchanged. Otherwise it is joined onto the
 * directory component of @p json_path, so the two files can be moved together as a directory.
 */
std::string resolve_sidx_path( std::string const& json_path, std::string const& sidx_filename );

#endif // SPEAR_COMMANDS_MAP_FORMAT_H_
