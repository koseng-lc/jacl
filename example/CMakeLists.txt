cmake_minimum_required(VERSION 2.8.3)
project(control_sim)

# QT
# Prevent Qt keywords with Python
ADD_DEFINITIONS(-DQT_NO_KEYWORDS)

set(Qt5Core_DIR /opt/qt59/lib/cmake/Qt5Core)
set(Qt5Gui_DIR /opt/qt59/lib/cmake/Qt5Gui)
set(Qt5Widgets_DIR /opt/qt59/lib/cmake/Qt5Widgets)

find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5Widgets REQUIRED)

# Main Program
set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)

add_executable(${PROJECT_NAME}
    # Header
    mainwindow.h
    controller_dialog.h
    # Source
    main.cpp
    mainwindow.cpp
    mainwindow.ui
    controller_dialog.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${Qt5Core_INCLUDE_DIRS}
    ${Qt5Gui_INCLUDE_DIRS}
    ${Qt5Widgets_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME}
    # Qt
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    # JACL
    jacl_lib
)
