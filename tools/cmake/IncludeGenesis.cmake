# SPEAR - Sorting Petabytes of Environmental and Ancient Reads
# Copyright (C) 2026 Lucas Czech
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact:
# Lucas Czech <lucas.czech@sund.ku.dk>
# University of Copenhagen, Globe Institute, Section for GeoGenetics
# Oster Voldgade 5-7, 1350 Copenhagen K, Denmark

# ------------------------------------------------------------------------------
#   Setup genesis
# ------------------------------------------------------------------------------

# Force to set the options used for genesis. Bit cumbersome in cmake...
# For now, spear is always compiled with C++20, so we deactivate
# the genesis auto detection of the highest available C++ standard here.
set(GENESIS_CPP_STANDARD "20" CACHE BOOL "Set genesis CPP standard." FORCE)
set(GENESIS_BUILD_APPLICATIONS OFF CACHE BOOL "" FORCE)
set(GENESIS_BUILD_PYTHON_MODULE OFF CACHE BOOL "" FORCE)
set(GENESIS_BUILD_TESTS OFF CACHE BOOL "" FORCE)

# Check if the submodule was initialized; otherwise fetch from GitHub
message (STATUS "Looking for genesis")
set(genesis_SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/deps/genesis)
if (EXISTS "${genesis_SUBMODULE_DIR}/CMakeLists.txt")
    message(STATUS "Using genesis submodule")
    set(genesis_SOURCE "${genesis_SUBMODULE_DIR}")
    add_subdirectory("${genesis_SUBMODULE_DIR}" "${CMAKE_BINARY_DIR}/genesis")
else()
    message(STATUS "Using genesis from GitHub")
    include(FetchContent)
    FetchContent_Declare(
        genesis
        GIT_REPOSITORY "https://github.com/lczech/genesis.git"
        GIT_TAG        "${genesis_GIT_TAG}"
    )
    FetchContent_MakeAvailable(genesis)
    set(genesis_SOURCE "${genesis_SOURCE_DIR}")
endif()

# Set the include dir for genesis headers
include_directories("${genesis_SOURCE}/lib")

# Now check targets exported by fetched version
if(TARGET genesis)
    message(STATUS "Fetched genesis with target 'genesis'")

elseif(TARGET genesis::genesis)
    message(STATUS "Fetched genesis with target 'genesis::genesis'")
    add_library(genesis ALIAS genesis::genesis)

else()
    # No exported target — fabricate a minimal one
    message(STATUS "No exported target found; creating INTERFACE genesis")
    add_library(genesis INTERFACE)
    target_include_directories(
        genesis INTERFACE
        $<BUILD_INTERFACE:${genesis_SOURCE}/lib>
        # $<INSTALL_INTERFACE:include>
    )
endif()

# Normalize to plain genesis
if(NOT TARGET genesis AND TARGET genesis::genesis)
    add_library(genesis ALIAS genesis::genesis)
endif()
if(TARGET genesis AND NOT TARGET genesis::genesis)
    add_library(genesis::genesis ALIAS genesis)
endif()

# Example usage:
# target_link_libraries(my_app PUBLIC genesis::genesis)
