#ifndef SPEAR_TOOLS_VERSION_H_
#define SPEAR_TOOLS_VERSION_H_

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

#include <string>

// =================================================================================================
//      SPEAR Version
// =================================================================================================

inline std::string spear_version()
{
    return "v0.0.1"; // #SPEAR_VERSION#
}

inline std::string spear_header()
{
    return "\
   .|'''.|  '||'''|. '||''''|     |     '||''|.          \n\
   ||..  '   ||   ||  ||  .      |||     ||   ||         \n\
    ''|||.   ||...|'  ||''|     |  ||    ||''|'          \n\
  .     '||  ||       ||       .''''|.   ||   |.         \n\
  '|....|'  .||.     .||....| .|.  .||. .||.  '|..       \n\
<<==================================================>>=＞\n\
   " + spear_version() + " (c) 2026 by Lucas Czech\n";
}

inline std::string spear_title()
{
    return "SPEAR: Sorting Petabytes of Environmental and Ancient Reads";
}

#endif // include guard
