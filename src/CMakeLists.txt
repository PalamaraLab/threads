# This file is part of https://github.com/PalamaraLab/threads which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Find boost
find_package(Boost REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

set(
        threads_infer_src
        Demography.cpp
        Node.cpp
        State.cpp
        Threads.cpp
        HMM.cpp
        TGEN.cpp
        ImputationMatcher.cpp
        Matcher.cpp
        ViterbiLowMem.cpp
        ThreadsLowMem.cpp
)

set(
        threads_infer_hdr
        Demography.hpp
        Node.hpp
        State.hpp
        Threads.hpp
        HMM.hpp
        TGEN.hpp
        TgenSegment.hpp
        ImputationMatcher.hpp
        Matcher.hpp
        ViterbiLowMem.hpp
        ThreadsLowMem.hpp
)

add_library(threads_infer STATIC ${threads_infer_src} ${threads_infer_hdr})

target_include_directories(
        threads_infer INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:>
)

set_target_properties(threads_infer PROPERTIES PUBLIC_HEADER "${threads_infer_hdr}")

target_link_libraries(threads_infer PRIVATE Boost::headers)
target_link_libraries(threads_infer PRIVATE project_warnings)
target_link_libraries(threads_infer PUBLIC Eigen3::Eigen)

if (PYTHON_BINDINGS)
    set_target_properties(threads_infer PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(threads_infer_python_bindings threads_infer_pybind.cpp)
    target_link_libraries(threads_infer_python_bindings PRIVATE threads_infer)
endif ()
