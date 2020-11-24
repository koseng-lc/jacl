cmake_minimum_required(VERSION 2.8.3)
project(basic_linalg)

add_executable(${PROJECT_NAME}1
    basic_linalg1.cpp
)

target_link_libraries(${PROJECT_NAME}1
    jacl_lib
)