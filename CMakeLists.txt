#
# JACL CMake
#
cmake_minimum_required(VERSION 2.8.3)
project(control_sim)

add_compile_options(-std=c++11)

# QT

# Prevent Qt keywords with Python
ADD_DEFINITIONS(-DQT_NO_KEYWORDS)

set(Qt5Core_DIR /opt/qt59/lib/cmake/Qt5Core)
set(Qt5Gui_DIR /opt/qt59/lib/cmake/Qt5Gui)
set(Qt5Widgets_DIR /opt/qt59/lib/cmake/Qt5Widgets)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)

find_package(PythonLibs 3.5 REQUIRED)

find_package(Boost 1.70.0 REQUIRED COMPONENTS
    system
    python
    numpy
    thread
    chrono
)

# Main Program

set(CMAKE_INCLUDE_CURRENT_DIR on)

include_directories(
        include
        src
        build

        ${PYTHON_INCLUDE_DIRS}

	${Qt5Core_INCLUDE_DIRS}
	${Qt5Gui_INCLUDE_DIRS}
	${Qt5Widgets_INCLUDE_DIRS}
)

# Python Binder

add_library(simulator SHARED include/simulator.h)

target_link_libraries(simulator
    Boost::system Boost::python
    ${PYTHON_LIBRARIES}
)
set_target_properties(simulator PROPERTIES SUFFIX .so)
set_target_properties(simulator PROPERTIES PREFIX "")
set_target_properties(simulator PROPERTIES LINKER_LANGUAGE CXX)

# JACL Lib

add_library(physical_parameter include/physical_parameter.h)
set_target_properties(physical_parameter PROPERTIES LINKER_LANGUAGE CXX)

add_library(state_space include/state_space.h)
set_target_properties(state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(state_space armadillo physical_parameter)

add_library(transfer_matrix include/transfer_matrix.h src/transfer_matrix.cpp)

add_library(observer include/observer.h src/observer.cpp)

add_library(simulator2 SHARED include/simulator2.h)

target_link_libraries(simulator2
    Boost::system Boost::python Boost::thread Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
)

set_target_properties(simulator2 PROPERTIES LINKER_LANGUAGE CXX)

add_library(jacl include/jacl)
set_target_properties(jacl PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(jacl state_space transfer_matrix simulator simulator2)

# Test

add_library(hello_ext SHARED src/hello.cpp)

target_link_libraries(hello_ext Boost::python)
set_target_properties(hello_ext PROPERTIES SUFFIX .so)
set_target_properties(hello_ext PROPERTIES PREFIX "")

add_executable(embed_test src/embed_test.cpp)
target_link_libraries(embed_test Boost::system Boost::python ${PYTHON_LIBRARIES})

# Main Program

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

add_executable(${PROJECT_NAME}
        # Header
	include/mainwindow.h
        # Source
	src/main.cpp
	src/mainwindow.cpp
	src/mainwindow.ui
)

target_link_libraries(${PROJECT_NAME}
        # Qt
	Qt5::Core
	Qt5::Gui
        Qt5::Widgets
        # JACL
        jacl
)
