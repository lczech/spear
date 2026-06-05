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

#include "spear/inverted_index/pfor.hpp"

#include <stdexcept>
#include <string>

extern "C" {
#include <turbopfor/ic.h>
}

namespace spear::inverted_index {

// =================================================================================================
//     32-bit
// =================================================================================================

std::size_t pfor_bound_32( std::size_t n ) noexcept
{
    return p4nbound32( n );
}

std::size_t pfor_encode_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nenc32( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4ndec32( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_delta_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4ndenc32( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_delta_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nddec32( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_delta1_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nd1enc32( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_delta1_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nd1dec32( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_zigzag_32(
    std::uint32_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nzenc32( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_zigzag_32(
    std::uint8_t* in,
    std::size_t n,
    std::uint32_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nzdec32( reinterpret_cast<unsigned char*>( in ), n, out );
}

// =================================================================================================
//     64-bit
// =================================================================================================

std::size_t pfor_bound_64( std::size_t n ) noexcept
{
    return p4nbound64( n );
}

std::size_t pfor_encode_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nenc64( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4ndec64( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_delta_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4ndenc64( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_delta_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nddec64( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_delta1_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nd1enc64( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_delta1_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nd1dec64( reinterpret_cast<unsigned char*>( in ), n, out );
}

std::size_t pfor_encode_zigzag_64(
    std::uint64_t* in,
    std::size_t n,
    std::uint8_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nzenc64( in, n, reinterpret_cast<unsigned char*>( out ));
}

std::size_t pfor_decode_zigzag_64(
    std::uint8_t* in,
    std::size_t n,
    std::uint64_t* out
) {
    if( n == 0 ) {
        return 0;
    }

    return p4nzdec64( reinterpret_cast<unsigned char*>( in ), n, out );
}

} // namespace spear::inverted_index
