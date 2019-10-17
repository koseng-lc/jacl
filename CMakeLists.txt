#
# JACL CMake
#
cmake_minimum_required(VERSION 2.8.3)
project(control_sim)

add_compile_options(-std=c++11)

set(Qt5Core_DIR /opt/qt59/lib/cmake/Qt5Core)
set(Qt5Gui_DIR /opt/qt59/lib/cmake/Qt5Gui)
set(Qt5Widgets_DIR /opt/qt59/lib/cmake/Qt5Widgets)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)

find_package(PythonLibs 3.5 REQUIRED)
find_package(PythonInterp 3.5 REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    python
    numpy
)

set(CMAKE_INCLUDE_CURRENT_DIR on)

include_directories(
        include
        src
        build
        /home/koseng/pylib/venv/include/python3.5m
	${Qt5Core_INCLUDE_DIRS}
	${Qt5Gui_INCLUDE_DIRS}
	${Qt5Widgets_INCLUDE_DIRS}
)
# Python Binder
add_library(hello_ext SHARED src/hello.cpp)

target_link_libraries(hello_ext boost_python-py35)
set_target_properties(hello_ext PROPERTIES SUFFIX .so)
set_target_properties(hello_ext PROPERTIES PREFIX "")

add_library(simulator SHARED include/simulator.h)

add_dependencies(simulator Boost)

target_link_libraries(simulator boost_python-py35 Boost::numpy ${PYTHON_LIBRARIES})
set_target_properties(simulator PROPERTIES SUFFIX .so)
set_target_properties(simulator PROPERTIES PREFIX "")
set_target_properties(simulator PROPERTIES LINKER_LANGUAGE CXX)

# JACL Lib
add_library(state_space include/state_space.h include/physical_parameter.h)
set_target_properties(state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(state_space armadillo)

add_library(transfer_matrix include/transfer_matrix.h src/transfer_matrix.cpp)

add_library(jacl include/jacl include/defs.h)
target_link_libraries(jacl state_space transfer_matrix)

# Main Program

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

add_executable(${PROJECT_NAME}
	include/mainwindow.h

	src/main.cpp
	src/mainwindow.cpp
	src/mainwindow.ui
)

target_link_libraries(${PROJECT_NAME}
	Qt5::Core
	Qt5::Gui
        Qt5::Widgets

        jacl
)
