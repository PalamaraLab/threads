# This file is part of https://github.com/PalamaraLab/TDPBWT which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Find boost
find_package(Boost) # COMPONENTS iostreams REQUIRED)

set(
        TDPBWT_src
        Demography.cpp
        Node.cpp
        State.cpp
        DPBWT.cpp
)

set(
        TDPBWT_hdr
        Demography.hpp
        Node.hpp
        State.hpp
        DPBWT.hpp
)

add_library(TDPBWT STATIC ${TDPBWT_src} ${TDPBWT_hdr})

target_include_directories(
        TDPBWT INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:>
)

set_target_properties(TDPBWT PROPERTIES PUBLIC_HEADER "${TDPBWT_hdr}")

# target_link_libraries(TDPBWT PRIVATE Boost::headers Boost::iostreams)
target_link_libraries(TDPBWT PRIVATE project_warnings)

# target_include_directories(TDPBWT PUBLIC ${Boost_INCLUDE_DIRS})

if (PYTHON_BINDINGS)
    set_target_properties(TDPBWT PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(TDPBWT_python_bindings TDPBWT_pybind.cpp)
    target_link_libraries(TDPBWT_python_bindings PRIVATE TDPBWT)
endif ()
