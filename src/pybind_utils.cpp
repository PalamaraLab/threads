// This file is part of the Threads software suite.
// Copyright (C) 2024-2025 Threads Developers.
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

#include "pybind_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

py::tuple threading_instructions_get_state(const ThreadingInstructions& ti) {
    // Instructions stored as a 4-tuple consisting of list of per-instruction
    // data, positions, start and end.
    py::list inst_list;
    for (const auto& inst : ti.instructions) {
        inst_list.append(
            py::make_tuple(
                inst.starts,
                inst.tmrcas,
                inst.targets,
                inst.mismatches
            )
        );
    }

    return py::make_tuple(inst_list, ti.positions, ti.start, ti.end);
}

ThreadingInstructions threading_instructions_set_state(py::tuple tup) {
    // Unpack the 4-tuple into threading instructions and positions
    std::vector<ThreadingInstruction> insts;
    py::list inst_list = tup[0];
    for (py::handle inst : inst_list) {
        // Each list entry must be a 4-tuple in same order as get_state above
        py::tuple inst_tup = inst.cast<py::tuple>();
        auto starts = inst_tup[0].cast<std::vector<int>>();
        auto tmrcas = inst_tup[1].cast<std::vector<double>>();
        auto targets = inst_tup[2].cast<std::vector<int>>();
        auto mismatches = inst_tup[3].cast<std::vector<int>>();
        insts.emplace_back(starts, tmrcas, targets, mismatches);
    }

    std::vector<int> positions = tup[1].cast<std::vector<int>>();
    ThreadingInstructions instructions{insts, positions};
    instructions.start = tup[2].cast<int>();
    instructions.end = tup[3].cast<int>();
    return instructions;
}
