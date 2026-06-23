# ------------------------------------------------------------------------------
#   Sources
# ------------------------------------------------------------------------------

# Collect library sources
file(GLOB_RECURSE spear_lib_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# ------------------------------------------------------------------------------
#   Build Object and Static Lib
# ------------------------------------------------------------------------------

# Object library: compiles the sources once into .o files
add_library(spear_lib_obj OBJECT ${spear_lib_sources})
target_compile_features(spear_lib_obj PUBLIC cxx_std_20)

target_include_directories(spear_lib_obj
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(spear_lib_obj
    PUBLIC
        genesis_lib_static
    PRIVATE
        turbopfor
        wfa2_static
)

# Static archive built from the object files
add_library(spear_lib_static STATIC
    $<TARGET_OBJECTS:spear_lib_obj>
)

target_compile_features(spear_lib_static PUBLIC cxx_std_20)
target_include_directories(spear_lib_static
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(spear_lib_static
    PUBLIC
        genesis_lib_static
        PRIVATE
        turbopfor
        wfa2_static
)
