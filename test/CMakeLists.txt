# This file is part of https://github.com/PalamaraLab/threads which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Get Catch2
Include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_SHALLOW TRUE
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(Catch2)
set(test_src
    test_demography.cpp
    test_viterbi_state.cpp
)

# Main unit test exe
add_executable(unit_tests
    ${test_src}
)

target_link_libraries(unit_tests
    PRIVATE
        threads_arg
        Catch2::Catch2WithMain
)

# automatically discover tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
