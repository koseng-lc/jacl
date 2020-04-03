cmake_minimum_required(VERSION 2.8.3)
project(example)

# QT
# Prevent Qt keywords with Python
ADD_DEFINITIONS(-DQT_NO_KEYWORDS)

set(Qt5Core_DIR /opt/qt59/lib/cmake/Qt5Core)
set(Qt5Gui_DIR /opt/qt59/lib/cmake/Qt5Gui)
set(Qt5Widgets_DIR /opt/qt59/lib/cmake/Qt5Widgets)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Boost REQUIRED COMPONENTS
    system
    python
    numpy
    thread
    chrono
)

# Main Program
set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

include_directories(${CMAKE_SOURCE_DIR}/include)

# Lib

# System Simulator
add_library(system_sim system_sim.h)
target_link_libraries(system_sim PRIVATE
    Boost::system Boost::python Boost::chrono
    ${PYTHON_LIBRARIES}
    jacl_lib
)
set_target_properties(system_sim PROPERTIES LINKER_LANGUAGE CXX)

# Controller
add_library(controller controller.h)
set_target_properties(controller PROPERTIES LINKER_LANGUAGE CXX)

# Controller Simulator
add_library(controller_sim controller_sim.h)
target_link_libraries(controller_sim PRIVATE jacl_lib)
set_target_properties(controller_sim PROPERTIES LINKER_LANGUAGE CXX)

# Position Controller Sim
add_library(posctrl_sim posctrl_sim.h)
target_link_libraries(posctrl_sim PRIVATE jacl_lib)
set_target_properties(posctrl_sim PROPERTIES LINKER_LANGUAGE CXX)

# Speed Controller Sim
add_library(spdctrl_sim spdctrl_sim.h)
target_link_libraries(spdctrl_sim PRIVATE jacl_lib)
set_target_properties(spdctrl_sim PROPERTIES LINKER_LANGUAGE CXX)

add_executable(control_sim
    # Header
    mainwindow.h
    controller_dialog.h
    # Source
    main.cpp
    mainwindow.cpp
    mainwindow.ui
    controller_dialog.cpp
)

target_include_directories(control_sim PRIVATE
    ${Qt5Core_INCLUDE_DIRS}
    ${Qt5Gui_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
)

target_link_libraries(control_sim
    # Qt
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    # JACL
    system_sim
    controller_sim
    posctrl_sim
    spdctrl_sim
)
