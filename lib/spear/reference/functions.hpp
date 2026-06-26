#ifndef SPEAR_LIB_REFERENCE_FUNCTIONS_H_
#define SPEAR_LIB_REFERENCE_FUNCTIONS_H_

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
 * @brief Free functions operating on ReferenceCollection: building from FASTA files and JSON I/O.
 *
 * @file
 * @ingroup reference
 */

#include "spear/reference/collection.hpp"
#include "genesis/util/format/json/document.hpp"

#include <string>
#include <vector>

namespace spear {
namespace reference {

// =================================================================================================
//     Builder
// =================================================================================================

/**
 * @brief Build a ReferenceCollection from a list of FASTA files.
 *
 * For each file, a FAI index (`.fai`) or sequence dictionary (`.dict`) is used if present,
 * avoiding a full sequence scan. A warning is logged if neither is available, since the fallback
 * reads the entire FASTA to obtain sequence names and lengths.
 *
 * If @p store_sites is true, the full sequence data is read from FASTA regardless of index
 * availability and stored in SequenceRecord::sites.
 *
 * Note: when loaded from FAI or dict, SequenceRecord::header equals the sequence name only
 * (the first whitespace-delimited token), since those formats do not preserve the full label.
 */
ReferenceCollection build_collection(
    std::vector<std::string> const& fasta_paths,
    bool store_sites = false
);

// =================================================================================================
//     JSON I/O
// =================================================================================================

/**
 * @brief Serialize a ReferenceCollection to a JSON document.
 *
 * Stores the FASTA file list and per-sequence header, length, and file_index.
 * Sequence `sites` data is omitted (too large for metadata); lengths are re-derived from the
 * prefix-sum offsets. The document can be round-tripped via collection_from_json().
 */
genesis::util::format::JsonDocument collection_to_json( ReferenceCollection const& collection );

/**
 * @brief Deserialize a ReferenceCollection from a JSON document produced by collection_to_json().
 *
 * Constructs and returns a new ReferenceCollection populated via the public add() interface.
 * Sequence `sites` data is not stored in JSON and will be empty in the returned collection.
 */
ReferenceCollection collection_from_json( genesis::util::format::JsonDocument const& doc );

} // namespace reference
} // namespace spear

#endif // SPEAR_LIB_REFERENCE_FUNCTIONS_H_
