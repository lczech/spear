#ifndef SPEAR_LIB_INVERTED_INDEX_STATS_H_
#define SPEAR_LIB_INVERTED_INDEX_STATS_H_

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
 * @brief
 *
 * @file
 * @ingroup util
 */

#include <cstdint>
#include <vector>

namespace spear::inverted_index {

// =================================================================================================
//     Inverted Index Stats
// =================================================================================================

/**
 * @brief Statistics about an Inverted Index.
 *
 * For instance produced by InvertedIndexBuilder::write().
 *
 * The posting_count_histogram has size max_postings_per_term + 2:
 *   - histogram[0]               : number of terms with no postings (empty)
 *   - histogram[1..max]          : number of terms with exactly that many postings
 *   - histogram[max + 1]         : number of terms that were capped
 *
 * The total_blob_bytes field reports the total number of bytes written for the compressed
 * posting lists (not counting the offset table or footer).
 */
struct InvertedIndexStats
{
    std::vector<std::uint64_t> posting_count_histogram;
    std::uint64_t total_blob_bytes = 0;
};

} // namespace spear::inverted_index

#endif // SPEAR_LIB_INVERTED_INDEX_STATS_H_
