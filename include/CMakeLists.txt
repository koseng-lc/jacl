#-- JACL Lib Project
cmake_minimum_required(VERSION 2.8.3)
project(jacl_lib)

add_compile_options(-std=c++14)

find_package(PythonLibs REQUIRED)
find_package(Boost REQUIRED COMPONENTS
    system
    python
    numpy
    thread
    chrono
)

set(JACL_HEADER_PREFIX ${PROJECT_SOURCE_DIR}/jacl)

include_directories(${PROJECT_SOURCE_DIR})

#-- Pattern
add_library(observer_pattern ${JACL_HEADER_PREFIX}/pattern/observer.hpp)
set_target_properties(observer_pattern PROPERTIES LINKER_LANGUAGE CXX)

#-- Physical Parameter
add_library(physical_parameter ${JACL_HEADER_PREFIX}/physical_parameter.hpp)
set_target_properties(physical_parameter PROPERTIES LINKER_LANGUAGE CXX)

#-- State Space
add_library(linear_state_space ${JACL_HEADER_PREFIX}/linear_state_space.hpp)
set_target_properties(linear_state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(linear_state_space PRIVATE armadillo physical_parameter observer_pattern)

add_library(nonlinear_state_space ${JACL_HEADER_PREFIX}/nonlinear_state_space.hpp)
set_target_properties(nonlinear_state_space PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(nonlinear_state_space PRIVATE armadillo physical_parameter observer_pattern)

#-- System
add_library(base_system ${JACL_HEADER_PREFIX}/system/base_system.hpp)
set_target_properties(base_system PROPERTIES LINKER_LANGUAGE CXX)
#target_include_directories(base_system PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(base_system PRIVATE armadillo physical_parameter)

#-- Continuous System
add_library(continuous_system ${JACL_HEADER_PREFIX}/system/continuous_system.hpp)
set_target_properties(continuous_system PROPERTIES LINKER_LANGUAGE CXX)
#target_include_directories(continuous_system PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(continuous_system PRIVATE base_system)

#-- Discrete System
add_library(discrete_system ${JACL_HEADER_PREFIX}/system/discrete_system.hpp)
set_target_properties(discrete_system PROPERTIES LINKER_LANGUAGE CXX)
#target_include_directories(discrete_system PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(discrete_system PRIVATE base_system)

#-- Observer
add_library(observer ${JACL_HEADER_PREFIX}/system/observer.hpp)
target_link_libraries(observer linear_state_space)
set_target_properties(observer PROPERTIES LINKER_LANGUAGE CXX)

add_library(continuous_observer ${JACL_HEADER_PREFIX}/system/continuous_observer.hpp)
target_link_libraries(continuous_observer linear_state_space)
set_target_properties(continuous_observer PROPERTIES LINKER_LANGUAGE CXX)

add_library(discrete_observer ${JACL_HEADER_PREFIX}/system/discrete_observer.hpp)
target_link_libraries(discrete_observer linear_state_space)
set_target_properties(discrete_observer PROPERTIES LINKER_LANGUAGE CXX)

#-- Transfer Matrix
add_library(transfer_matrix ${JACL_HEADER_PREFIX}/transfer_matrix.hpp)
set_target_properties(transfer_matrix PROPERTIES LINKER_LANGUAGE CXX)

#-- Traits
add_library(traits ${JACL_HEADER_PREFIX}/traits.hpp)
set_target_properties(traits PROPERTIES LINKER_LANGUAGE CXX)

#-- Linear Algebra
add_library(linear_algebra ${JACL_HEADER_PREFIX}/linear_algebra.hpp)
target_link_libraries(linear_algebra PRIVATE armadillo)
set_target_properties(linear_algebra PROPERTIES LINKER_LANGUAGE CXX)

#-- Numerical Methods
add_library(numerical_methods ${JACL_HEADER_PREFIX}/numerical_methods.hpp)
set_target_properties(numerical_methods PROPERTIES LINKER_LANGUAGE CXX)

#-- LTI Common
add_library(lti_common ${JACL_HEADER_PREFIX}/lti_common.hpp)
target_link_libraries(lti_common PRIVATE armadillo)
set_target_properties(lti_common PROPERTIES LINKER_LANGUAGE CXX)

#-- Simulator(more precisely is plotter)
add_library(sim ${JACL_HEADER_PREFIX}/sim.hpp)
target_link_libraries(sim PRIVATE
    Boost::system Boost::python Boost::thread Boost::chrono
    ${PYTHON_LIBRARIES}
    linear_state_space
)
set_target_properties(sim PROPERTIES LINKER_LANGUAGE CXX)

#-- Simulator(more precisely is plotter)
add_library(plotter ${JACL_HEADER_PREFIX}/plotter.hpp)
target_link_libraries(plotter PRIVATE
    Boost::system Boost::python Boost::thread Boost::chrono
    ${PYTHON_LIBRARIES}
    continuous_system
    discrete_system
)
set_target_properties(plotter PROPERTIES LINKER_LANGUAGE CXX)

#-- LFT
add_library(lft ${JACL_HEADER_PREFIX}/lft.hpp)
target_link_libraries(lft PRIVATE linear_state_space)
set_target_properties(lft PROPERTIES LINKER_LANGUAGE CXX)

#-- Upper LFT
add_library(upper_lft ${JACL_HEADER_PREFIX}/upper_lft.hpp)
target_link_libraries(upper_lft PRIVATE linear_state_space lft)
set_target_properties(upper_lft PROPERTIES LINKER_LANGUAGE CXX)

#-- Lower LFT
add_library(lower_lft ${JACL_HEADER_PREFIX}/lower_lft.hpp)
target_link_libraries(lower_lft PRIVATE linear_state_space lft)
set_target_properties(lower_lft PROPERTIES LINKER_LANGUAGE CXX)

#-- Pole Placement
add_library(pole_placement ${JACL_HEADER_PREFIX}/pole_placement.hpp)
target_link_libraries(pole_placement PRIVATE Boost::system Boost::thread Boost::chrono)
set_target_properties(pole_placement PROPERTIES LINKER_LANGUAGE CXX)

#-- ARE
add_library(are ${JACL_HEADER_PREFIX}/are.hpp)
target_link_libraries(are PRIVATE Boost::system Boost::python ${Python_LIBRARIES} linear_state_space)
set_target_properties(are PROPERTIES LINKER_LANGUAGE CXX)

#-- H-infinity Synthesis
add_library(h_inf ${JACL_HEADER_PREFIX}/synthesis/h_inf.hpp)
target_link_libraries(h_inf PRIVATE upper_lft are lti_common numerical_methods)
set_target_properties(h_inf PROPERTIES LINKER_LANGUAGE CXX)

#-- Instrument Fault Detection
add_library(ifd ${JACL_HEADER_PREFIX}/ifd.hpp)
target_link_libraries(ifd PRIVATE continuous_system discrete_system)
set_target_properties(ifd PROPERTIES LINKER_LANGUAGE CXX)

#-- Parser
add_library(parser ${JACL_HEADER_PREFIX}/parser.hpp)
set_target_properties(parser PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(parser PRIVATE linear_state_space traits)

set(JACL_LIB_LIST
    physical_parameter
    linear_state_space
    nonlinear_state_space
    transfer_matrix
    traits
    linear_algebra
    numerical_methods
    lti_common
    sim
    plotter
    lft
    upper_lft
    lower_lft
    pole_placement
    are
    h_inf
    observer_pattern
    parser
    base_system
    continuous_system
    discrete_system
    observer
    continuous_observer
    discrete_observer
    ifd
)

#-- JACL Libraries
add_library(jacl_lib ${JACL_HEADER_PREFIX}/jacl.hpp)
set_target_properties(jacl_lib PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(jacl_lib
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${JACL_INCLUDE_PATH}>
    $<INSTALL_INTERFACE:${PYTHON_INCLUDE_DIRS}>
)
target_link_libraries(jacl_lib PRIVATE ${JACL_LIB_LIST})

#-- Install
install(TARGETS jacl_lib ${JACL_LIB_LIST} EXPORT jacl_lib_export
        ARCHIVE DESTINATION ${JACL_LIB_PATH}
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/jacl/ DESTINATION ${JACL_HEADER_PATH}/)
