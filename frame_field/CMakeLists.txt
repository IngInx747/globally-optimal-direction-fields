project(libffgen)

file(GLOB SRCS
    "*.h"
    "*.hh"
    "*.hpp"
    "*.c"
    "*.cc"
    "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

# mesh
target_link_libraries(${PROJECT_NAME} libmesh)
add_dependencies(${PROJECT_NAME} libmesh)

# numerical methods
target_link_libraries(${PROJECT_NAME} libnm)
add_dependencies(${PROJECT_NAME} libnm)

# export
target_include_directories(${PROJECT_NAME} PUBLIC .)
