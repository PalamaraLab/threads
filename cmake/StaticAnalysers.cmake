# This file is part of https://github.com/PalamaraLab/arg_needle_lib which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Based on https://github.com/lefticus/cpp_starter_project

option(ENABLE_CPPCHECK "Enable static analysis with cppcheck" OFF)
option(ENABLE_CLANG_TIDY "Enable static analysis with clang-tidy" OFF)

if(ENABLE_CPPCHECK)
  find_program(CPPCHECK cppcheck)
  if(CPPCHECK)
    set(CMAKE_CXX_CPPCHECK
        ${CPPCHECK}
        --suppress=missingInclude
        --suppress=unusedFunction
        --suppress=unmatchedSuppression
        --enable=all
        --inconclusive)
    message(STATUS "Using cppcheck ${CPPCHECK}")
  else()
    message(SEND_ERROR "cppcheck requested but executable not found")
  endif()
endif()

if(ENABLE_CLANG_TIDY)
  find_program(CLANGTIDY NAMES clang-tidy-11 clang-tidy-10 clang-tidy-9 clang-tidy-8 clang-tidy )
  if(CLANGTIDY)
    set(CMAKE_CXX_CLANG_TIDY ${CLANGTIDY} -extra-arg=-Wno-unknown-warning-option)
    message(STATUS "Using clang-tidy ${CLANGTIDY}")
  else()
    message(SEND_ERROR "clang-tidy requested but executable not found")
  endif()
endif()
