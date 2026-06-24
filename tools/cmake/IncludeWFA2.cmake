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
#   Make available WFA2
# ------------------------------------------------------------------------------

# Check if the submodule was initialized; otherwise fetch from GitHub.
# We build from sources directly to avoid their CMakeLists.txt requiring pkg-config.
message(STATUS "Looking for WFA2")
set(WFA2_SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/deps/wfa2)
if(EXISTS "${WFA2_SUBMODULE_DIR}/wavefront/wavefront_align.h")
    message(STATUS "Using WFA2 submodule")
    set(WFA2_SOURCE     "${WFA2_SUBMODULE_DIR}")
    set(WFA2_PARENT_DIR "${PROJECT_SOURCE_DIR}/deps")
else()
    message(STATUS "Using WFA2 from GitHub")
    FetchContent_Declare(
        wfa2
        GIT_REPOSITORY "https://github.com/smarco/WFA2-lib.git"
        GIT_TAG        "${wfa2_GIT_TAG}"
        SOURCE_DIR     "${CMAKE_BINARY_DIR}/deps/wfa2"
    )
    # Use Populate (not MakeAvailable) to download without processing WFA2's own
    # CMakeLists.txt, which defines wfa2_static and would collide with ours below.
    FetchContent_GetProperties(wfa2)
    if(NOT wfa2_POPULATED)
        FetchContent_Populate(wfa2)
    endif()
    set(WFA2_SOURCE     "${wfa2_SOURCE_DIR}")
    set(WFA2_PARENT_DIR "${CMAKE_BINARY_DIR}/deps")
endif()

# ------------------------------------------------------------------------------
#   Setup WFA2
# ------------------------------------------------------------------------------

if(NOT EXISTS "${WFA2_SOURCE}/wavefront/wavefront_align.h")
    message(FATAL_ERROR
        "WFA2 not found at ${WFA2_SOURCE}. "
        "Did you initialize the submodule?"
    )
endif()

# Glob all C sources from the relevant upstream subdirectories.
file(GLOB wfa2_SOURCES
    "${WFA2_SOURCE}/wavefront/*.c"
    "${WFA2_SOURCE}/alignment/*.c"
    "${WFA2_SOURCE}/system/*.c"
    "${WFA2_SOURCE}/utils/*.c"
)

add_library(wfa2_static STATIC ${wfa2_SOURCES})
add_library(wfa2::wfa2 ALIAS wfa2_static)

target_compile_features(wfa2_static PUBLIC c_std_99)
set_target_properties(wfa2_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

# External includes: consumers use #include "wfa2/..." via the parent dir.
# WFA2 root also PUBLIC so WFA2's own headers can resolve their internal cross-includes
# (e.g. wfa.h includes "system/mm_allocator.h" relative to the wfa2 root).
target_include_directories(wfa2_static SYSTEM PUBLIC
    "${WFA2_PARENT_DIR}"
    "${WFA2_SOURCE}"
)

# Suppress warnings from this third-party library
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(wfa2_static PRIVATE -w)
endif()

target_compile_definitions(wfa2_static PRIVATE NDEBUG)
