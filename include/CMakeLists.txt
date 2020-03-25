#-- JACL Lib
cmake_minimum_required(VERSION 2.8.3)
project(jacl_lib)

add_compile_options(-std=c++14)

find_package(PythonLibs REQUIRED)

find_package(Boost REQUIRED COMPONENTS
    system
    python3
    numpy
    thread
    chrono
)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PYTHON_INCLUDE_DIRS}
)

set(JACL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(JACL_HEADER_PREFIX ${JACL_SOURCE_DIR}/jacl)

# Physical Parameter
add_library(physical_parameter ${JACL_HEADER_PREFIX}/physical_parameter.h)
set_target_properties(physical_parameter PROPERTIES LINKER_LANGUAGE CXX)

# State Space
add_library(state_space ${JACL_HEADER_PREFIX}/state_space.h)
set_target_properties(state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(state_space armadillo physical_parameter)

# Transfer Matrix
add_library(transfer_matrix ${JACL_HEADER_PREFIX}/transfer_matrix.h ${JACL_HEADER_PREFIX}/transfer_matrix.cpp)

# Traits
add_library(traits ${JACL_HEADER_PREFIX}/traits.h)
set_target_properties(traits PROPERTIES LINKER_LANGUAGE CXX)

# Linear Algebra
add_library(linear_algebra ${JACL_HEADER_PREFIX}/linear_algebra.h)
target_link_libraries(linear_algebra armadillo)
set_target_properties(linear_algebra PROPERTIES LINKER_LANGUAGE CXX)

# Numerical Methods
add_library(numerical_methods ${JACL_HEADER_PREFIX}/numerical_methods.h)
set_target_properties(numerical_methods PROPERTIES LINKER_LANGUAGE CXX)

# LTI Common
add_library(lti_common ${JACL_HEADER_PREFIX}/lti_common.h)
target_link_libraries(lti_common armadillo)
set_target_properties(lti_common PROPERTIES LINKER_LANGUAGE CXX)

# Observer
add_library(observer ${JACL_HEADER_PREFIX}/observer.h)
set_target_properties(observer PROPERTIES LINKER_LANGUAGE CXX)

# Simulator(more precisely is plotter)
add_library(sim SHARED ${JACL_HEADER_PREFIX}/sim.h)
target_link_libraries(sim
    Boost::system Boost::python3 Boost::thread Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
)
set_target_properties(sim PROPERTIES LINKER_LANGUAGE CXX)

# Observer Simulator
add_library(observer_sim ${JACL_HEADER_PREFIX}/observer_sim.h)
target_link_libraries(observer_sim
    Boost::system Boost::python3 Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
    sim
    observer
)
set_target_properties(observer_sim PROPERTIES LINKER_LANGUAGE CXX)

# LFT
add_library(lft ${JACL_HEADER_PREFIX}/lft.h)
target_link_libraries(lft state_space)
set_target_properties(lft PROPERTIES LINKER_LANGUAGE CXX)

# Upper LFT
add_library(upper_lft ${JACL_HEADER_PREFIX}/upper_lft.h)
target_link_libraries(upper_lft PUBLIC state_space lft)
set_target_properties(upper_lft PROPERTIES LINKER_LANGUAGE CXX)

# Lower LFT
add_library(lower_lft ${JACL_HEADER_PREFIX}/lower_lft.h)
target_link_libraries(lower_lft PUBLIC state_space lft)
set_target_properties(lower_lft PROPERTIES LINKER_LANGUAGE CXX)

# Pole Placement
add_library(pole_placement ${JACL_HEADER_PREFIX}/pole_placement.h)
target_link_libraries(pole_placement Boost::system Boost::thread Boost::chrono)
set_target_properties(pole_placement PROPERTIES LINKER_LANGUAGE CXX)

# ARE
add_library(are ${JACL_HEADER_PREFIX}/are.h)
target_link_libraries(are Boost::system Boost::python3 ${Python_LIBRARIES} state_space)
set_target_properties(are PROPERTIES LINKER_LANGUAGE CXX)

# H-infinity Synthesis
add_library(h_infinity_synthesis ${JACL_HEADER_PREFIX}/h_infinity_synthesis.h)
target_link_libraries(h_infinity_synthesis upper_lft are lti_common numerical_methods)
set_target_properties(h_infinity_synthesis PROPERTIES LINKER_LANGUAGE CXX)

# System Simulator
add_library(system_sim ${JACL_HEADER_PREFIX}/system_sim.h)
target_link_libraries(system_sim
    Boost::system Boost::python3 Boost::chrono
    ${PYTHON_LIBRARIES}
    state_space
    sim
)
set_target_properties(system_sim PROPERTIES LINKER_LANGUAGE CXX)

# Controller
add_library(controller ${JACL_HEADER_PREFIX}/controller.h)
set_target_properties(controller PROPERTIES LINKER_LANGUAGE CXX)

# Controller Simulator
add_library(controller_sim ${JACL_HEADER_PREFIX}/controller_sim.h)
target_link_libraries(controller_sim
    sim
)
set_target_properties(controller_sim PROPERTIES LINKER_LANGUAGE CXX)

# Position Controller Sim
add_library(posctrl_sim ${JACL_HEADER_PREFIX}/posctrl_sim.h)
target_link_libraries(posctrl_sim
    sim
)
set_target_properties(posctrl_sim PROPERTIES LINKER_LANGUAGE CXX)

# Speed Controller Sim
add_library(spdctrl_sim ${JACL_HEADER_PREFIX}/spdctrl_sim.h)
target_link_libraries(spdctrl_sim
    sim
)
set_target_properties(spdctrl_sim PROPERTIES LINKER_LANGUAGE CXX)

# JACL
add_library(jacl_lib ${JACL_HEADER_PREFIX}/jacl.hpp)
set_target_properties(jacl_lib PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(jacl_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(jacl_lib traits state_space transfer_matrix sim h_infinity_synthesis controller_sim)
