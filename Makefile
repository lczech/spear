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
#     This makefile wraps around CMake in order to make the use straight forward.
#     A simple call to "make" suffices to build the whole of SPEAR.
#
#     Build type is controlled by which target you call first (or explicitly):
#       make          — RELEASE by default on a fresh build; preserves existing config otherwise
#       make release  — (re)configure as RELEASE, then compile
#       make debug    — (re)configure as DEBUG (assertions, debug symbols), then compile
#
#     This script is mainly intended for fast development, as it 'misuses' CMake
#     directly as a build system instead of a build system _generator_.
# --------------------------------------------------------------------------------------------------

# Run CMake if not yet done or if CMakeLists.txt has changed, then compile.
# Does not pass -DCMAKE_BUILD_TYPE so that an existing cache is not overridden;
# on a fresh build CMakeLists.txt defaults to RELEASE when the variable is unset.
all: build/CMakeCache.txt
	@echo "Running make..."
	$(MAKE) -s -C build
	@echo "Done."
.PHONY: all

# Run CMake if not yet done or if CMakeLists.txt has changed.
build/CMakeCache.txt: CMakeLists.txt
	@echo "Running CMake..."
	@mkdir -p build
	@cd build && cmake ..

# Explicitly (re)configure as RELEASE and compile.
release:
	@echo "Running CMake (release)..."
	@mkdir -p build
	@cd build && cmake .. -DCMAKE_BUILD_TYPE=RELEASE
	@echo "Running make..."
	$(MAKE) -s -C build
	@echo "Done."
.PHONY: release

# Explicitly (re)configure as DEBUG (assertions, debug symbols) and compile.
debug:
	@echo "Running CMake (debug)..."
	@mkdir -p build
	@cd build && cmake .. -DCMAKE_BUILD_TYPE=DEBUG
	@echo "Running make..."
	$(MAKE) -s -C build
	@echo "Done."
.PHONY: debug

# Special make that also includes new files.
# We first touch all inner CMake files so that their glob search for files is rerun.
# This ensures that all new files are compiled, even when doing incremental builds.
update: build/CMakeCache.txt
	@echo "Running make with new files..."
	@touch CMakeLists.txt
	@if [ -d libs/genesis ]; then touch libs/genesis/CMakeLists.txt; fi
	$(MAKE) -s -C build
.PHONY: update

# Clean up all build targets.
clean:
	@echo "Running clean..."
	@rm -rf bin
	@rm -rf build
.PHONY: clean
