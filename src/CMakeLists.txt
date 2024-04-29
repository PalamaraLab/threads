# This file is part of https://github.com/PalamaraLab/threads which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Find boost
find_package(Boost)
find_package(Eigen3)

set(
        threads_src
        Demography.cpp
        Node.cpp
        State.cpp
        ThreadsFastLS.cpp
        HMM.cpp
        TGEN.cpp
        ImputationMatcher.cpp
        Matcher.cpp
        ViterbiLowMem.cpp
        ThreadsLowMem.cpp
)

set(
        threads_hdr
        Demography.hpp
        Node.hpp
        State.hpp
        ThreadsFastLS.hpp
        HMM.hpp
        TGEN.hpp
        ImputationMatcher.hpp
        Matcher.hpp
        ViterbiLowMem.hpp
        ThreadsLowMem.hpp
)

add_library(threads STATIC ${threads_src} ${threads_hdr})

target_include_directories(
        threads INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:>
)

set_target_properties(threads PROPERTIES PUBLIC_HEADER "${threads_hdr}")

target_link_libraries(threads PRIVATE Boost::headers)
target_link_libraries(threads PUBLIC Eigen3::Eigen)
target_link_libraries(threads PRIVATE project_warnings)

if (PYTHON_BINDINGS)
    set_target_properties(threads PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(threads_python_bindings threads_pybind.cpp)
    target_link_libraries(threads_python_bindings PRIVATE threads)
endif ()
