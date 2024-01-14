cmake_minimum_required(VERSION 3.15)
project(argweaver)

# Require C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set install prefix to project root
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

# I prefer #ifndef NDEBUG, but argweaver uses #ifdef DEBUG
add_compile_definitions($<$<CONFIG:Debug>:DEBUG>)

# Source files
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
file(GLOB SOURCES "${SOURCE_DIR}/argweaver/*.cpp")

# Compile argweaver as a static library
add_library(obj OBJECT ${SOURCES})
include_directories(${SOURCE_DIR})
add_library(argweaver STATIC $<TARGET_OBJECTS:obj>)

# All executables are in src/
file(GLOB EXECUTABLE_SOURCES "${SOURCE_DIR}/*.cpp")

# Compile executables
foreach(executable_cpp ${EXECUTABLE_SOURCES})
    get_filename_component(executable ${executable_cpp} NAME_WE)
    add_executable(${executable} ${executable_cpp})
    target_link_libraries(${executable} argweaver)
    install(TARGETS ${executable})
endforeach()
