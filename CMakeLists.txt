#
# JACL CMake
#

cmake_minimum_required(VERSION 2.8.3)
project(jacl)

set(JACL_PACKAGE_VERSION_MAJOR 1)
set(JACL_PACKAGE_VERSION_MINOR 0)
set(JACL_PACKAGE_VERSION_PATCH 0)

set(CMAKE_INCLUDE_CURRENT_DIR on)

add_subdirectory(example)
add_subdirectory(include)
