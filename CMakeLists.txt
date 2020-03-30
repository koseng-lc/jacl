#
# JACL CMake
#

cmake_minimum_required(VERSION 2.8.3)
project(jacl)

set(CMAKE_BUILD_TYPE Release)

set(JACL_PACKAGE_VERSION_MAJOR 1)
set(JACL_PACKAGE_VERSION_MINOR 0)
set(JACL_PACKAGE_VERSION_PATCH 0)
set(JACL_PACKAGE_VERSION "${JACL_PACKAGE_VERSION_MAJOR}.${JACL_PACKAGE_VERSION_MINOR}.${JACL_PACKAGE_VERSION_PATCH}")

set(JACL_INSTALL_PREFIX "/media/koseng/4A7AE1C07AE1A941/lit/capstone/capstone_ws") # Change as you need
set(JACL_MAIN_INSTALL_PATH "${JACL_INSTALL_PREFIX}/jacl-${JACL_PACKAGE_VERSION}")
set(JACL_INCLUDE_PATH "${JACL_MAIN_INSTALL_PATH}/include")
set(JACL_HEADER_PATH "${JACL_INCLUDE_PATH}/jacl")
set(JACL_LIB_PATH "${JACL_MAIN_INSTALL_PATH}/lib/jacl")
set(JACL_CMAKE_CONFIG_INSTALL_PATH "${JACL_MAIN_INSTALL_PATH}/share/jacl/cmake")

# Fix this !!!
find_package(PythonLibs REQUIRED)

add_subdirectory(example)
add_subdirectory(include)

# for build-tree
set(CONF_INCLUDE_DIRS "${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}")
set(CONF_TARGETS_DIR "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTOY}/Export/_${JACL_CMAKE_CONFIG_INSTALL_PATH}")
configure_file(jacl-config.cmake.in
    "${CMAKE_BINARY_DIR}/jacl-config.cmake" @ONLY)
# for install-tree
set(CONF_INCLUDE_DIRS "${JACL_INCLUDE_PATH}")
set(CONF_TARGETS_DIR "${JACL_CMAKE_CONFIG_INSTALL_PATH}")
configure_file(jacl-config.cmake.in
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/jacl-config.cmake" @ONLY)

install(FILES "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/jacl-config.cmake" DESTINATION ${JACL_CMAKE_CONFIG_INSTALL_PATH})
install(EXPORT jacl_lib_export DESTINATION ${JACL_CMAKE_CONFIG_INSTALL_PATH}
        FILE jacl-targets.cmake
)
