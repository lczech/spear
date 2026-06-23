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
#   Make available TurboPFor
# ------------------------------------------------------------------------------

# Check if the submodule was initialized; otherwise fetch from GitHub
message (STATUS "Looking for TurboPFor")
set(TURBOPFOR_SUBMODULE_DIR ${PROJECT_SOURCE_DIR}/deps/TurboPFor)
if (EXISTS "${TURBOPFOR_SUBMODULE_DIR}/include/ic.h")
    message(STATUS "Using TurboPFor submodule")
    set(TURBOPFOR_SOURCE "${TURBOPFOR_SUBMODULE_DIR}")
else()
    message(STATUS "Using TurboPFor from GitHub")
    FetchContent_Declare(
        TurboPFor
        GIT_REPOSITORY "https://github.com/powturbo/TurboPFor-Integer-Compression"
        GIT_TAG        "${TurboPFor_GIT_TAG}"
        SOURCE_DIR     "${CMAKE_BINARY_DIR}/TurboPFor"
    )
    FetchContent_MakeAvailable(TurboPFor)
    set(TURBOPFOR_SOURCE "${turbopfor_SOURCE_DIR}")
endif()

# ------------------------------------------------------------------------------
#   Setup TurboPFor
# ------------------------------------------------------------------------------

if(NOT EXISTS "${TURBOPFOR_SOURCE}/include/ic.h")
    message(FATAL_ERROR
        "TurboPFor not found at ${TURBOPFOR_SOURCE}. "
        "Did you initialize the submodule?"
    )
endif()

# Base library sources, matching upstream LIB=...
set(TURBOPFOR_BASE_SOURCES
    "${TURBOPFOR_SOURCE}/lib/bic.c"
    "${TURBOPFOR_SOURCE}/lib/bitunpack.c"
    "${TURBOPFOR_SOURCE}/lib/bitpack.c"
    "${TURBOPFOR_SOURCE}/lib/bitutil.c"
    "${TURBOPFOR_SOURCE}/lib/eliasfano.c"
    "${TURBOPFOR_SOURCE}/lib/fp.c"
    "${TURBOPFOR_SOURCE}/lib/transpose.c"
    "${TURBOPFOR_SOURCE}/lib/transpose_.c"
    "${TURBOPFOR_SOURCE}/lib/trlec.c"
    "${TURBOPFOR_SOURCE}/lib/trled.c"
    "${TURBOPFOR_SOURCE}/lib/vp4c.c"
    "${TURBOPFOR_SOURCE}/lib/vp4d.c"
    "${TURBOPFOR_SOURCE}/lib/v8.c"
    "${TURBOPFOR_SOURCE}/lib/v8pack.c"
    "${TURBOPFOR_SOURCE}/lib/vint.c"
    "${TURBOPFOR_SOURCE}/lib/vsimple.c"
    "${TURBOPFOR_SOURCE}/lib/vbit.c"
)

add_library(turbopfor STATIC ${TURBOPFOR_BASE_SOURCES})
add_library(turbopfor::turbopfor ALIAS turbopfor)
# set_source_files_properties(${TURBOPFOR_BASE_SOURCES}
#     PROPERTIES COMPILE_OPTIONS "-mno-avx2;-mno-bmi;-mno-bmi2"
# )

target_compile_features(turbopfor PUBLIC c_std_99)
set_target_properties(turbopfor PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Public headers plus wrapper include directory for <turbopfor/ic.h>
set(TURBOPFOR_WRAPPER_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${TURBOPFOR_WRAPPER_DIR}/turbopfor")
file(GENERATE OUTPUT "${TURBOPFOR_WRAPPER_DIR}/turbopfor/ic.h"
    CONTENT "#pragma once\n#include <ic.h>\n"
)

target_include_directories(turbopfor
    PUBLIC
    "${TURBOPFOR_SOURCE}/include"
    "${TURBOPFOR_WRAPPER_DIR}"
)

# Upstream uses these broadly.
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(turbopfor PRIVATE
        -fstrict-aliasing
        -w
        -Wall
        -pedantic
    )
endif()

# GCC 14+ made -Wincompatible-pointer-types a hard error for C; suppress it
# for this third-party library without touching its sources.
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_compile_options(turbopfor PRIVATE -Wno-incompatible-pointer-types)
endif()

target_compile_definitions(turbopfor PRIVATE NDEBUG)

# Match upstream arch defaults reasonably closely.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(turbopfor PRIVATE
            -mssse3
            -mno-avx2
            -mno-bmi
            -mno-bmi2
        )
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)$")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(turbopfor PRIVATE -march=armv8-a)
    endif()
    target_compile_definitions(turbopfor PRIVATE _NAVX2)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(ppc64le)$")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(turbopfor PRIVATE -mcpu=power9 -mtune=power9)
    endif()
    target_compile_definitions(turbopfor PRIVATE __SSSE3__)
    target_compile_definitions(turbopfor PRIVATE _NAVX2)
else()
    target_compile_definitions(turbopfor PRIVATE _NAVX2)
endif()

# Linux upstream links m and rt.
if(UNIX AND NOT APPLE)
    target_link_libraries(turbopfor PUBLIC m rt)
elseif(APPLE)
    target_link_libraries(turbopfor PUBLIC m)
endif()

# --------------------------------------------------------------------------
# Extra AVX2 objects on x86_64
# --------------------------------------------------------------------------

# vp4c.c vp4d.c transpose.c bitpack.c bitunpack.c bitutil.c
# Important:
# - compile the SAME sources a second time with AVX2 enabled
# - the sources themselves expose distinct ...256v... symbols

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
    message (STATUS "TurboPFor with AVX2 support on x86_64")

    set(TURBOPFOR_AVX2_SOURCES
        "${TURBOPFOR_SOURCE}/lib/vp4c.c"
        "${TURBOPFOR_SOURCE}/lib/vp4d.c"
        "${TURBOPFOR_SOURCE}/lib/transpose.c"
        "${TURBOPFOR_SOURCE}/lib/bitpack.c"
        "${TURBOPFOR_SOURCE}/lib/bitunpack.c"
        "${TURBOPFOR_SOURCE}/lib/bitutil.c"
    )

    add_library(turbopfor_avx2 OBJECT ${TURBOPFOR_AVX2_SOURCES})

    target_compile_features(turbopfor_avx2 PRIVATE c_std_99)
    set_target_properties(turbopfor_avx2 PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_include_directories(turbopfor_avx2 PRIVATE
        "${TURBOPFOR_SOURCE}/include"
    )

    target_compile_definitions(turbopfor_avx2 PRIVATE NDEBUG)

    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(turbopfor_avx2 PRIVATE
            -O3
            -fstrict-aliasing
            -w
            -Wall
            -pedantic
            -mavx2
            -mbmi
            -mbmi2
        )

        # Force AVX2/BMI/BMI2 on these specific cloned sources.
        # set_source_files_properties(${TURBOPFOR_AVX2_SOURCES} PROPERTIES
        #     COMPILE_OPTIONS "-mavx2;-mbmi;-mbmi2"
        # )
    endif()

    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        target_compile_options(turbopfor_avx2 PRIVATE
            -Wno-incompatible-pointer-types
        )
    endif()

    target_sources(turbopfor PRIVATE $<TARGET_OBJECTS:turbopfor_avx2>)
endif()
