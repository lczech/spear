deps
====

External dependencies of spear. We include them here as submodules, such that they can be downloaded once via `git init --recursive`, for instance for deploying spear on systems without internet access. However, our CMake setup also fetches these internally when they are not yet present, for ease of use.
