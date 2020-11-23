cmake_minimum_required(VERSION 2.8.3)
project(basic_linalg)

add_executable(${PROJECT_NAME}_1
    basic_linalg1.cpp
)

target_link_libraries(${PROJECT_NAME}_1
    jacl_lib
)