# CMakeLists of libraries to build as static libs.
# The libs are used by multiple executables and other libs.

# OpenMesh
set(OPENMESH_DOCS false BOOL)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/openmesh")

# Eigen
add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE "${CMAKE_CURRENT_LIST_DIR}/eigen")
