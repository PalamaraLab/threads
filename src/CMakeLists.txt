# This file is part of https://github.com/PalamaraLab/threads which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Find Boost
find_package(Boost REQUIRED)
message(STATUS "Found Boost ${Boost_VERSION}")

# When using GCC C++17 and Boost ICL < 1.65, need to turn off new-ttp-matching.
# See https://github.com/PalamaraLab/TDPBWT/issues/7
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(Boost_FOUND AND Boost_VERSION VERSION_LESS "1.65.0")
        message(STATUS "GCC is used with Boost version less than 1.65, disabling -fnew-ttp-matching.")
        add_compile_options(-fno-new-ttp-matching)
    endif()
endif()

# Try to find Eigen3 version 3.4 locally installed
find_package(Eigen3 3.4 QUIET)

if(NOT Eigen3_FOUND)
    include(FetchContent)
    message(STATUS "Suitable local version of Eigen3 (>= 3.4) not found, fetching from repository...")
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        9df21dc8b4b576a7aa5c0094daa8d7e8b8be60f0  # 3.4 branch at 2024-03-15
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
    set(Eigen3_VERSION "3.4")
endif()
message(STATUS "Found Eigen3 ${Eigen3_VERSION}")

# Threads static library
set(threads_arg_src
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

set(threads_arg_hdr
    Demography.hpp
    Node.hpp
    State.hpp
    ThreadsFastLS.hpp
    HMM.hpp
    TGEN.hpp
    TgenSegment.hpp
    ImputationMatcher.hpp
    Matcher.hpp
    ViterbiLowMem.hpp
    ThreadsLowMem.hpp
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
    PUBLIC
        Eigen3::Eigen
)

# Conditionally create python bindings
if(PYTHON_BINDINGS)
    set_target_properties(threads_arg
        PROPERTIES
            POSITION_INDEPENDENT_CODE TRUE
    )

    pybind11_add_module(threads_arg_python_bindings
        threads_arg_pybind.cpp
    )

    target_link_libraries(threads_arg_python_bindings
        PRIVATE
            threads_arg
            project_warnings
    )
endif()
