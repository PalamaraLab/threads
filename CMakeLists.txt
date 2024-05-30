# This file is part of https://github.com/PalamaraLab/threads which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")

project(threads_arg LANGUAGES CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project settings including default build type
include(cmake/ProjectSettings.cmake)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitiser options if supported by compiler
include(cmake/Sanitisers.cmake)
enable_sanitisers(project_settings)

# allow for static analysis options
include(cmake/StaticAnalysers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)

option(PYTHON_BINDINGS "Whether to build the python bindings" ON)
if(PYTHON_BINDINGS)
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG aa304c9c7d725ffb9d10af08a3b34cb372307020 # Version 2.10.0
    )
    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
endif()

add_subdirectory(src)

option(ENABLE_TESTING "Enable test builds" OFF)
if(ENABLE_TESTING)
    message(STATUS "Test builds enabled")
    enable_testing()
    add_subdirectory(test)
endif()
