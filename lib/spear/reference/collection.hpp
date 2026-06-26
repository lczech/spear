#ifndef SPEAR_LIB_REFERENCE_COLLECTION_H_
#define SPEAR_LIB_REFERENCE_COLLECTION_H_

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
 * @brief Ordered collection of reference sequences with a shared global coordinate space.
 *
 * @file
 * @ingroup reference
 */

#include "genesis/sequence/sequence.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace spear {
namespace reference {

// =================================================================================================
//     ReferenceCollection
// =================================================================================================

/**
 * @brief Ordered collection of reference sequences spanning a shared global coordinate space.
 *
 * Sequences are added incrementally (possibly from multiple FASTA files) and form a concatenated
 * coordinate space: sequence i occupies global positions
 * [global_offset_of(i), global_offset_of(i) + sequences()[i].length).
 * This matches the coordinate convention used by InvertedIndexBuilder, where genome bin positions
 * are offsets into this concatenated space.
 *
 * Internally an n+1 prefix-sum sentinel vector is maintained: offsets_[0] = 0,
 * offsets_[i+1] = offsets_[i] + length_i, so offsets_.back() == total_length().
 *
 * Use build_collection() (see functions.hpp) to populate from FASTA files with automatic
 * FAI / sequence-dictionary detection. Use to_json() / from_json() (see functions.hpp) to
 * serialize and deserialize the collection.
 */
class ReferenceCollection
{
public:

    // -------------------------------------------------------------------------
    //     SequenceRecord
    // -------------------------------------------------------------------------

    /**
     * @brief Metadata and optional sequence data for a single reference sequence entry.
     *
     * `header` stores the full FASTA label line (everything after `>`).
     * When loaded from a FAI or dict file, `header` equals the sequence name only,
     * since those formats do not preserve the full FASTA description.
     * `sites` is non-empty only when the collection was built with `store_sites = true`.
     */
    struct SequenceRecord
    {
        // Full FASTA label line (everything after '>').
        std::string header;

        // Index into fasta_files() identifying the source input file.
        std::size_t file_index = 0;

        // Length in bases of this sequence.
        std::size_t length = 0;

        // Actual nucleotide sequence; empty unless the collection was built with store_sites=true.
        std::string sites;

        /**
         * @brief Return the first whitespace-delimited token of `header`.
         *
         * This is the canonical sequence name as used by alignment tools and index formats.
         */
        std::string name() const
        {
            auto const end = header.find_first_of( " \t" );
            return ( end == std::string::npos ) ? header : header.substr( 0, end );
        }
    };

    // -------------------------------------------------------------------------
    //     FindResult
    // -------------------------------------------------------------------------

    /**
     * @brief Result returned by find_by_global_position().
     *
     * `sequence` points into the collection's internal vector and is valid for the
     * ReferenceCollection's lifetime. `local_offset` is the distance from the sequence's start
     * to the queried position (i.e., global_pos - global_offset_of(index)).
     */
    struct FindResult
    {
        // Pointer into sequences(); valid for the ReferenceCollection's lifetime.
        SequenceRecord const* sequence = nullptr;

        // Index of the sequence within sequences().
        std::size_t index = 0;

        // Offset of the queried position within the sequence: global_pos - sequence_start.
        std::size_t local_offset = 0;
    };

    // -------------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------------

    ReferenceCollection()  = default;
    ~ReferenceCollection() = default;

    ReferenceCollection( ReferenceCollection const& ) = default;
    ReferenceCollection( ReferenceCollection&& )      = default;

    ReferenceCollection& operator=( ReferenceCollection const& ) = default;
    ReferenceCollection& operator=( ReferenceCollection&& )      = default;

    // -------------------------------------------------------------------------
    //     Modifiers
    // -------------------------------------------------------------------------

    /**
     * @brief Add a sequence given its full FASTA label and length; no sequence data is stored.
     */
    void add( std::string const& fasta_path, std::string const& header, std::size_t length );

    /**
     * @brief Add a sequence given its full FASTA label and nucleotide bases; @p sites is moved in.
     *
     * The length is derived from `sites.size()`.
     */
    void add( std::string const& fasta_path, std::string const& header, std::string sites );

    /**
     * @brief Add a sequence from a genesis Sequence object, copying label and sites.
     */
    void add( std::string const& fasta_path, genesis::sequence::Sequence const& sequence );

    /**
     * @brief Add a sequence from a genesis Sequence object, moving the sites string.
     *
     * Avoids copying the (potentially large) sequence data when the caller no longer needs it.
     */
    void add( std::string const& fasta_path, genesis::sequence::Sequence&& sequence );

    /**
     * @brief Remove all sequences and file registrations, resetting to an empty collection.
     */
    void clear();

    // -------------------------------------------------------------------------
    //     Accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Return true if no sequences have been added.
     */
    bool empty() const {
        return sequences_.empty();
    }

    /**
     * @brief Return the number of sequences in the collection.
     */
    std::size_t size() const {
        return sequences_.size();
    }

    /**
     * @brief Return the SequenceRecord at @p index without bounds checking.
     */
    SequenceRecord const& operator[]( std::size_t index ) const {
        return sequences_[index];
    }

    /**
     * @brief Return the SequenceRecord at @p index with bounds checking.
     */
    SequenceRecord const& at( std::size_t index ) const {
        return sequences_.at( index );
    }

    /**
     * @brief Return the list of source FASTA file paths, in order of first appearance.
     */
    std::vector<std::string> const& fasta_files() const {
        return fasta_files_;
    }

    /**
     * @brief Return the source FASTA file path at @p index.
     */
    std::string const& fasta_file( std::size_t index ) const {
        return fasta_files_.at( index );
    }

    /**
     * @brief Return all SequenceRecord entries.
     */
    std::vector<SequenceRecord> const& sequences() const {
        return sequences_;
    }

    /**
     * @brief Return the global start offset of the sequence at @p index.
     *
     * The sequence occupies global positions
     * [global_offset_of(i), global_offset_of(i) + sequences()[i].length).
     */
    std::size_t global_offset_of( std::size_t index ) const {
        return offsets_[index];
    }

    /**
     * @brief Return the total number of bases across all sequences.
     */
    std::size_t total_length() const {
        return offsets_.back();
    }

    // -------------------------------------------------------------------------
    //     Position Lookup
    // -------------------------------------------------------------------------

    /**
     * @brief Return the sequence containing @p global_pos and the local offset within it.
     *
     * Uses binary search (O(log n)) on the prefix-sum offset table.
     *
     * @throws std::out_of_range if the collection is empty or global_pos >= total_length().
     */
    FindResult find_by_global_position( std::size_t global_pos ) const;

    // -------------------------------------------------------------------------
    //     Iterators
    // -------------------------------------------------------------------------

    using const_iterator = std::vector<SequenceRecord>::const_iterator;

    /**
     * @brief Return an iterator to the first SequenceRecord.
     */
    const_iterator begin() const { return sequences_.cbegin(); }

    /**
     * @brief Return a past-the-end iterator over SequenceRecord entries.
     */
    const_iterator end()   const { return sequences_.cend(); }

private:

    // Register fasta_path in fasta_files_ if it differs from the last seen file; return its index.
    std::size_t register_file_( std::string const& fasta_path );

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

    std::vector<std::string>    fasta_files_;
    std::vector<SequenceRecord> sequences_;

    // n+1 prefix-sum sentinel: offsets_[0]=0, offsets_[i+1]=offsets_[i]+length_i.
    // Initialized with {0} so total_length()==0 on an empty collection.
    std::vector<std::size_t>  offsets_ = { 0 };
};

} // namespace reference
} // namespace spear

#endif // SPEAR_LIB_REFERENCE_COLLECTION_H_
