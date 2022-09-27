# This file is part of https://github.com/PalamaraLab/arg_needle_lib which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")

project(threads LANGUAGES CXX VERSION 0.1)

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

# set(arg_needle_testdata_dir ${CMAKE_CURRENT_SOURCE_DIR}/testdata)
# add_definitions(-DARG_NEEDLE_TESTDATA_DIR=\"${arg_needle_testdata_dir}\")

option(PYTHON_BINDINGS "Whether to build the python bindings" ON)
if (PYTHON_BINDINGS)
    include(FetchContent)
    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG ffa346860b306c9bbfb341aed9c14c067751feb8 # Version 2.9.1
    )
    FetchContent_GetProperties(pybind11)
    if (NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif ()
endif ()

add_subdirectory(src)
# add_subdirectory(example)

option(ENABLE_TESTING "Enable Test Builds" OFF)
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()
