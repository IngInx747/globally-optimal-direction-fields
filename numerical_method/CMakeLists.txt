project(libnm)

file(GLOB SRCS
    "*.h"
    "*.hh"
    "*.hpp"
    "*.c"
    "*.cc"
    "*.cpp")

add_library(${PROJECT_NAME} ${SRCS})

# Eigen
target_link_libraries(${PROJECT_NAME} eigen)
add_dependencies(${PROJECT_NAME} eigen)

# export
target_include_directories(${PROJECT_NAME} PUBLIC .)
