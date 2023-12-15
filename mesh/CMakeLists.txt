project(libmesh)

file(GLOB SRCS
    "*.h"
    "*.hh"
    "*.hpp"
    "*.c"
    "*.cc"
    "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

# OpenMesh
target_link_libraries(${PROJECT_NAME} OpenMeshCore)
add_dependencies(${PROJECT_NAME} OpenMeshCore)

# export
target_include_directories(${PROJECT_NAME} PUBLIC .)
