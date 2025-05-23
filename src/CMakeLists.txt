# This file is part of the Threads software suite.
# Copyright (C) 2024-2025 Threads Developers.
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

# Find Boost
find_package(Boost REQUIRED)
message(STATUS "Found Boost ${Boost_VERSION}")

# Threads static library
set(threads_arg_src
    Demography.cpp
    Node.cpp
    State.cpp
    ThreadsFastLS.cpp
    HMM.cpp
    ImputationMatcher.cpp
    Matcher.cpp
    ViterbiLowMem.cpp
    ThreadsLowMem.cpp
    ThreadingInstructions.cpp
    DataConsistency.cpp
    AlleleAges.cpp
    GenotypeIterator.cpp
    VCFWriter.cpp
)

set(threads_arg_hdr
    Demography.hpp
    Node.hpp
    State.hpp
    ThreadsFastLS.hpp
    HMM.hpp
    TgenSegment.hpp
    ImputationMatcher.hpp
    Matcher.hpp
    ViterbiLowMem.hpp
    ThreadsLowMem.hpp
    ThreadingInstructions.hpp
    DataConsistency.hpp
    AlleleAges.hpp
    GenotypeIterator.hpp
    VCFWriter.hpp
)

add_library(threads_arg STATIC
    ${threads_arg_src}
)

target_include_directories(threads_arg INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:>
)

set_target_properties(threads_arg
    PROPERTIES
        PUBLIC_HEADER "${threads_arg_hdr}"
)

target_link_libraries(threads_arg
    PRIVATE
        Boost::headers
        project_warnings
)

# Conditionally create python bindings
if(PYTHON_BINDINGS)
    set_target_properties(threads_arg
        PROPERTIES
            POSITION_INDEPENDENT_CODE TRUE
    )

    pybind11_add_module(threads_arg_python_bindings
        threads_arg_pybind.cpp
        pybind_utils.cpp
    )

    target_link_libraries(threads_arg_python_bindings
        PRIVATE
            threads_arg
            project_warnings
    )

    install(TARGETS threads_arg_python_bindings LIBRARY DESTINATION .)
endif()
