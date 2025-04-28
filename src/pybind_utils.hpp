// This file is part of the Threads software suite.
// Copyright (C) 2025 Threads Developers.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef THREADS_ARG_PYBIND_UTILS_HPP
#define THREADS_ARG_PYBIND_UTILS_HPP

#include "ThreadingInstructions.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// get_state and set_state are used to pickle ThreadingInstructions so they
// may be moved between processes.
pybind11::tuple threading_instructions_get_state(const ThreadingInstructions& ti);
ThreadingInstructions threading_instructions_set_state(pybind11::tuple tup);

#endif // THREADS_ARG_PYBIND_UTILS_HPP
