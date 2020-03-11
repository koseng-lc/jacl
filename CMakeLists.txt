#
# JACL CMake
#

cmake_minimum_required(VERSION 2.8.3)
project(control_sim)

add_compile_options(-std=c++14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#message("${CMAKE_CXX_FLAGS}")

# QT
# Prevent Qt keywords with Python
ADD_DEFINITIONS(-DQT_NO_KEYWORDS)

set(Qt5Core_DIR /opt/qt59/lib/cmake/Qt5Core)
set(Qt5Gui_DIR /opt/qt59/lib/cmake/Qt5Gui)
set(Qt5Widgets_DIR /opt/qt59/lib/cmake/Qt5Widgets)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)

find_package(PythonLibs REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    system
    python3
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
    Boost::system Boost::python3
    ${PYTHON_LIBRARIES}
)
set_target_properties(simulator PROPERTIES SUFFIX .so)
set_target_properties(simulator PROPERTIES PREFIX "")
set_target_properties(simulator PROPERTIES LINKER_LANGUAGE CXX)

#-- JACL Lib

# Physical Parameter
add_library(physical_parameter include/physical_parameter.h)
set_target_properties(physical_parameter PROPERTIES LINKER_LANGUAGE CXX)

# State Space
add_library(state_space include/state_space.h)
set_target_properties(state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(state_space armadillo physical_parameter)

# Transfer Matrix
add_library(transfer_matrix include/transfer_matrix.h src/transfer_matrix.cpp)

# Traits
add_library(traits include/traits.h)
set_target_properties(traits PROPERTIES LINKER_LANGUAGE CXX)

# Linear Algebra
add_library(linear_algebra include/linear_algebra.h)
target_link_libraries(linear_algebra armadillo)
set_target_properties(linear_algebra PROPERTIES LINKER_LANGUAGE CXX)

# LTI Common
add_library(lti_common include/lti_common.h)
target_link_libraries(lti_common armadillo)
set_target_properties(lti_common PROPERTIES LINKER_LANGUAGE CXX)

# Observer
add_library(observer include/observer.h)
set_target_properties(observer PROPERTIES LINKER_LANGUAGE CXX)

# Simulator(more precisely is plotter)
add_library(sim SHARED include/sim.h)
target_link_libraries(sim
    Boost::system Boost::python3 Boost::thread Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
#    X11
)
set_target_properties(sim PROPERTIES LINKER_LANGUAGE CXX)

# Observer Simulator
add_library(observer_sim include/observer_sim.h)
target_link_libraries(observer_sim
    Boost::system Boost::python3 Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
    sim
    observer
)
set_target_properties(observer_sim PROPERTIES LINKER_LANGUAGE CXX)

# LFT
add_library(lft include/lft.h)
target_link_libraries(lft state_space)
set_target_properties(lft PROPERTIES LINKER_LANGUAGE CXX)

# Upper LFT
add_library(upper_lft include/upper_lft.h)
target_link_libraries(upper_lft PUBLIC state_space lft)
set_target_properties(upper_lft PROPERTIES LINKER_LANGUAGE CXX)

# Lower LFT
add_library(lower_lft include/lower_lft.h)
target_link_libraries(lower_lft PUBLIC state_space lft)
set_target_properties(lower_lft PROPERTIES LINKER_LANGUAGE CXX)

# Pole Placement
add_library(pole_placement include/pole_placement.h)
target_link_libraries(pole_placement Boost::system Boost::thread Boost::chrono)
set_target_properties(pole_placement PROPERTIES LINKER_LANGUAGE CXX)

# ARE
add_library(are include/are.h)
target_link_libraries(are state_space)
set_target_properties(are PROPERTIES LINKER_LANGUAGE CXX)

# H-infinity Synthesis
add_library(h_infinity_synthesis include/h_infinity_synthesis.h)
target_link_libraries(h_infinity_synthesis upper_lft are lti_common)
set_target_properties(h_infinity_synthesis PROPERTIES LINKER_LANGUAGE CXX)

# System Simulator
add_library(system_sim include/system_sim.h)
target_link_libraries(system_sim
    Boost::system Boost::python3 Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
    sim
)
set_target_properties(system_sim PROPERTIES LINKER_LANGUAGE CXX)

# Controller
add_library(controller include/controller.h)
set_target_properties(controller PROPERTIES LINKER_LANGUAGE CXX)

# Controller Simulator
add_library(controller_sim include/controller_sim.h)
target_link_libraries(controller_sim
    Boost::system Boost::python3 Boost::chrono
    ${Python_LIBRARIES}
    sim
)
set_target_properties(controller_sim PROPERTIES LINKER_LANGUAGE CXX)

# JACL
add_library(jacl include/jacl)
set_target_properties(jacl PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(jacl traits state_space transfer_matrix simulator sim h_infinity_synthesis controller_sim)

#-- End of JACL Lib part

# Test

add_library(hello_ext SHARED src/hello.cpp)

target_link_libraries(hello_ext Boost::python3)
set_target_properties(hello_ext PROPERTIES SUFFIX .so)
set_target_properties(hello_ext PROPERTIES PREFIX "")

add_executable(embed_test src/embed_test.cpp)
target_link_libraries(embed_test Boost::system Boost::python3 ${PYTHON_LIBRARIES})

# Main Program

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

add_executable(${PROJECT_NAME}
        # Header
	include/mainwindow.h
        include/controller_dialog.h
        # Source
	src/main.cpp
	src/mainwindow.cpp
	src/mainwindow.ui
        src/controller_dialog.cpp
)

target_link_libraries(${PROJECT_NAME}
        # Qt
	Qt5::Core
	Qt5::Gui
        Qt5::Widgets
        # JACL
        jacl
)
