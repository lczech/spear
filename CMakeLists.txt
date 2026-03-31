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

# --------------------------------------------------------------------------------------------------
#   CMake Init
# --------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)

# Fun with colors!
if (NOT WIN32)
    string(ASCII 27 Esc)

    set(ColorBlue   "${Esc}[34m") # Build information.
    set(ColorGreen  "${Esc}[32m") # Optional additions.
    set(ColorYellow "${Esc}[33m") # Tipps for the user.
    set(ColorRed    "${Esc}[31m") # Warnings and Errors.
    set(ColorEnd    "${Esc}[0m")
endif()

set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set (CMAKE_DISABLE_SOURCE_CHANGES  ON)

# Need to make sure that AppleClang is actually recognized and not mistakingly identifyed
# as normal Clang. See https://cmake.org/cmake/help/v3.0/policy/CMP0025.html
# This is important for the test cmake script, where we need to set some warnings for clang,
# but not for apple clang.
if( NOT ( CMAKE_VERSION VERSION_LESS 3.0 ))
    # Appleclang compatibility
    cmake_policy(SET CMP0025 NEW)
    # Support for :: aliasing
    cmake_policy(SET CMP0028 NEW)
endif()

# --------------------------------------------------------------------------------------------------
#   Project Definition and Requirements
# --------------------------------------------------------------------------------------------------

project ( spear CXX )

# We already check here whether the requirements for autotools are met,
# as we need it in spear, but want to avoid the misleading error message from genesis
# in case we don't have autotools.
find_program(AUTOCONF_RUNTIME autoconf)
if( AUTOCONF_RUNTIME STREQUAL "AUTOCONF_RUNTIME-NOTFOUND" )

    message (
        STATUS "${ColorRed}We could not find autotools/autoconf, "
        "which is needed to compile htslib.${ColorEnd}"
    )
    IF(APPLE)
        message ( STATUS "On MacOS, try: `brew install autoconf automake libtool`" )
    ELSE()
        message ( STATUS "On Linux, try: `sudo apt-get install autoconf autotools-dev`" )
    ENDIF()
    message( FATAL_ERROR "Required autotools/autoconf for building htslib not found.")
endif()

# Some fun for the users to figure this out :-( Hope to improve this error message in the future.
if(GENESIS_USE_HTSLIB STREQUAL "OFF")
    message(
        FATAL_ERROR "${ColorRed}Support for htslib has been deactivated, probably because some "
        "previous attempt to build spear recommended to set `-DGENESIS_USE_HTSLIB=OFF`. "
        "Unfortunately, this recommendation is produced by our library genesis, which can work "
        "without htslib, and hence produces this message. For spear however, we do need htslib. "
        "Hence, please instead install the dependencies required to build htslib, and try again. "
        "Sorry for the confusion.${ColorEnd}"
    )
endif()

# --------------------------------------------------------------------------------------------------
#   Build Options
# --------------------------------------------------------------------------------------------------

# The build type is determined by the standard cmake variable CMAKE_BUILD_TYPE.
# Thus, if you want to change the build type, call `cmake .. -DCMAKE_BUILD_TYPE=DEBUG`
IF(CMAKE_BUILD_TYPE STREQUAL "")
    set (CMAKE_BUILD_TYPE RELEASE)
ENDIF()

# Print build type. Useful for debugging user issues.
message (STATUS "${ColorBlue}spear build type: ${CMAKE_BUILD_TYPE}${ColorEnd}")

# Store everything in `bin`.
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )

# Option to enable building a fully static version, that is, statical linking of system libraries.
# This works with GCC, but might cause trouble with Clang, so it is disabled by default.
# Most users who compile on their own probaby won't need the fully portable version.
# Genesis on the other hand is always linked statically, in order to be able to move around the
# final binary at least on the same system.
option ( SPEAR_BUILD_STATIC  "Build the static version of the executable." OFF )

if( SPEAR_BUILD_STATIC AND ( NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" ))
    string(ASCII 27 Esc)
    message( STATUS "${Esc}[31m"
        "Compiler is not GCC. This might cause trouble with statically linking system libraries "
        "on Linux systems. If the compilation fails with some linking error, try a different compiler, "
        "or try to not build statically by deactivating the CMake option 'SPEAR_BUILD_STATIC'. "
    "${Esc}[0m" )
endif()

# --------------------------------------------------------------------------------------------------
#   Compiler and Linker Options
# --------------------------------------------------------------------------------------------------

if (MSVC OR (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC"))
    # MinGW/MSVC do not support much. Also, the above string comparison in the condition does
    # not work there, hence the extra `MSVC` check...
    set (WARN_FLAGS "-Wall")
else ()
    # For other operating systems, set the warn flags to a very high level - except unknown pragmas,
    # as this causes needless warnings with some pragma-based techniques.
    set (WARN_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -pedantic -pedantic-errors")

    # Furthermore, in debug mode, we also want all warnings to be errors, so that we do not miss them
    # In release, this might hinder users to successfully build genesis if any warnings remain
    # or occur only on their compiler, so let's avoid trouble for them.
    IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
        set (WARN_FLAGS "${WARN_FLAGS} -Werror")
    ENDIF()
endif ()

set (CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++20 ${WARN_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG   "-O2 -DDEBUG -g -ggdb3 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" )

IF(SPEAR_BUILD_STATIC)
    # Use different ways of telling CMake to build static, just to be sure.
    SET(BUILD_SHARED_LIBRARIES OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
    link_libraries("-static")
ENDIF()
message( STATUS "Static linking of system libraries: ${SPEAR_BUILD_STATIC}")

# Policy CMP0042: Enable MACOSX_RPATH
if (APPLE)
   set(CMAKE_MACOSX_RPATH ON)
endif (APPLE)

# Genesis uses LTO, and so do we.
if( NOT ( CMAKE_VERSION VERSION_LESS 3.9 ))
    cmake_policy(SET CMP0069 NEW)
endif()

# --------------------------------------------------------------------------------------------------
#   Download Dependencies
# --------------------------------------------------------------------------------------------------

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/tools/cmake" )

# Downloads dependencies if they are not there. We use a special script for this.
# include( "${CMAKE_CURRENT_LIST_DIR}/tools/cmake/DownloadDependency.cmake" )

# Store the commit hashes of our dependencies.
# These are replaced by tools/cmake/update_dependencies.sh to the hashes that are currently checked out.
# Thus, do not replace the hashes manually!
SET( CLI11_GIT_TAG "37bb6edc5317e99af72ef48405e65d9ca5218861" ) #CLI11_GIT_TAG#
SET( genesis_GIT_TAG "dadd370f13a5333611e5c078ce9c86ecee41cd9f" ) #genesis_GIT_TAG#

# Include the libaries and make them available here
include( IncludeCLI11 )
include( IncludeGenesis )

# --------------------------------------------------------------------------------------------------
#   Dependencies Settings
# --------------------------------------------------------------------------------------------------

# Add Genesis as dependency.
# add_subdirectory(${PROJECT_SOURCE_DIR}/libs/genesis)

# Use everything that Genesis exports, just to be sure that we use the same setup.
# add_definitions( ${GENESIS_DEFINITIONS} )
# include_directories( ${GENESIS_INCLUDE_DIR} )
# set( CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}          ${GENESIS_C_FLAGS}")
# set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${GENESIS_CXX_FLAGS}" )
# set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )

# Add CLI11 as dependency. It's header-only, so this is easy.
include_directories( SYSTEM "libs/CLI11/include" )

# ------------------------------------------------------------------------------
#   Sources
# ------------------------------------------------------------------------------

# Find all source files. The CMake documentation does not recommend this,
# but for our purposes, this is good enough for now.
file (GLOB_RECURSE spear_sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Add our own code. The main is only one cpp file, the rest are headers included from there.
# include_directories( "src" )
add_executable( spear ${spear_sources} )
target_include_directories(spear PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(spear PRIVATE genesis_lib_static)

# Link it against Genesis, and against all dependencies of Genesis.
# target_link_libraries( spear ${GENESIS_LINK_LIBRARIES} )
# target_link_libraries( spear PUBLIC genesis::genesis )

# We have to use LTO for spear as well for some compilers to properly link, I think.
if(GENESIS_HAS_LTO)
    message(STATUS "Building with LTO/IPO support")
    set_property(TARGET spear PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "Building without LTO/IPO support")
endif()

message( STATUS "CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}" )
message( STATUS "GENESIS_LINK_LIBRARIES ${GENESIS_LINK_LIBRARIES}" )
