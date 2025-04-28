# This file is part of the Threads software suite.
# Copyright (C) 2025 Threads Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")

project(threads_arg LANGUAGES CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_NEWPYTHON ON)

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
        GIT_TAG a2e59f0e7065404b44dfe92a28aca47ba1378dc4 # Version 2.13.6
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
