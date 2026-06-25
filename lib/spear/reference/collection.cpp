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

#include "spear/reference/collection.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace spear {
namespace reference {

// =================================================================================================
//     Modifiers
// =================================================================================================

void ReferenceCollection::add(
    std::string const& fasta_path,
    std::string const& header,
    std::uint64_t length
) {
    if( length == 0 ) {
        throw std::invalid_argument(
            "ReferenceCollection::add: sequence '" + header + "' has length zero"
        );
    }
    auto const fi = register_file_( fasta_path );
    sequences_.push_back({ header, fi, length, {} });
    offsets_.push_back( offsets_.back() + length );
}

void ReferenceCollection::add(
    std::string const& fasta_path,
    std::string const& header,
    std::string sites
) {
    if( sites.empty() ) {
        throw std::invalid_argument(
            "ReferenceCollection::add: sequence '" + header + "' has empty sites"
        );
    }
    auto const length = static_cast<std::uint64_t>( sites.size() );
    auto const fi = register_file_( fasta_path );
    sequences_.push_back({ header, fi, length, std::move( sites ) });
    offsets_.push_back( offsets_.back() + length );
}

void ReferenceCollection::add(
    std::string const& fasta_path,
    genesis::sequence::Sequence const& sequence
) {
    add( fasta_path, sequence.label(), sequence.sites() );
}

void ReferenceCollection::add(
    std::string const& fasta_path,
    genesis::sequence::Sequence&& sequence
) {
    add( fasta_path, sequence.label(), std::move( sequence.sites() ) );
}

void ReferenceCollection::clear()
{
    fasta_files_.clear();
    sequences_.clear();
    offsets_ = { 0 };
}

// =================================================================================================
//     Position Lookup
// =================================================================================================

ReferenceCollection::FindResult ReferenceCollection::find_by_global_position(
    std::uint64_t global_pos
) const {
    if( sequences_.empty() ) {
        throw std::out_of_range(
            "ReferenceCollection::find_by_global_position: collection is empty"
        );
    }
    if( global_pos >= total_length() ) {
        throw std::out_of_range(
            "ReferenceCollection::find_by_global_position: position " +
            std::to_string( global_pos ) + " >= total_length " +
            std::to_string( total_length() )
        );
    }

    // upper_bound on the prefix-sum table finds the first offset strictly greater than global_pos;
    // stepping back one gives the sequence whose range [offsets_[i], offsets_[i+1]) contains it.
    // offsets_[0] == 0, so the result is always at least begin()+1, making the step-back safe.
    auto const it = std::upper_bound( offsets_.begin(), offsets_.end(), global_pos );
    auto const idx = static_cast<std::size_t>( std::distance( offsets_.begin(), it ) ) - 1;

    return FindResult{ &sequences_[idx], idx, global_pos - offsets_[idx] };
}

// =================================================================================================
//     Private Helpers
// =================================================================================================

std::size_t ReferenceCollection::register_file_( std::string const& fasta_path )
{
    if( fasta_files_.empty() || fasta_files_.back() != fasta_path ) {
        fasta_files_.push_back( fasta_path );
    }
    return fasta_files_.size() - 1;
}

} // namespace reference
} // namespace spear
