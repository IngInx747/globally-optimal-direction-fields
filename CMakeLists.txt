# cmake_minimum_require must be the first command of the file
cmake_minimum_required(VERSION 3.0)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(FFGEN VERSION 1.0 LANGUAGES CXX)

if (WIN32)
    add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
endif ()

# we add the sub-directories that we want CMake to scan
add_subdirectory(3rdparty)
add_subdirectory(mesh)
add_subdirectory(numeric)
add_subdirectory(frame_field)
add_subdirectory(test)
