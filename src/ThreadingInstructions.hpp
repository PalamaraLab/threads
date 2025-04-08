// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
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

#ifndef THREADS_ARG_THREADING_INSTRUCTIONS_HPP
#define THREADS_ARG_THREADING_INSTRUCTIONS_HPP

#include "ViterbiLowMem.hpp"

#include <limits>
#include <sstream>
#include <vector>


class ThreadingInstruction {
// This is a container class for the instructions inferred for a single sequence
public:
    // Segment starts
    std::vector<int> starts;
    // Coalescence times (times to the most recent common ancestor)
    std::vector<double> tmrcas;
    // Closest cousins
    std::vector<int> targets;
    // Mismatches with closest cousin
    std::vector<int> mismatches;
    std::size_t num_segments = 0;
    std::size_t num_mismatches = 0;
    ThreadingInstruction(
        const std::vector<int>& _starts,
        const std::vector<double>& _tmrcas,
        const std::vector<int>& _targets,
        const std::vector<int>& _mismatches
    );
};

class ThreadingInstructionIterator {
// This is a container for iterating through a ThreadingInstruction object site-by-site
private:
    const ThreadingInstruction& instruction;
    const std::vector<int>& positions;
    int sites_processed = 0;
    int next_segment_start = 0;
    int next_mismatch = 0;
public:
    int num_sites = 0;
    int current_mismatch = 0;
    int current_segment = 0;
    int current_target = 0;
    double current_tmrca = 0;
    bool is_mismatch = false;
public:
    void increment_site(const int pos);
    ThreadingInstructionIterator(const ThreadingInstruction& _instruction, const std::vector<int>& _positions);
};

class ThreadingInstructions {
// This is a container class for threading instructions inferred by Threads
public:
    // Getters
    std::vector<std::vector<int>> all_starts();
    std::vector<std::vector<double>> all_tmrcas();
    std::vector<std::vector<int>> all_targets();
    std::vector<std::vector<int>> all_mismatches();

    // Constructors
    ThreadingInstructions(const std::vector<ThreadingInstruction>& _instructions, const std::vector<int>& _positions);
    ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int start, const int end, const std::vector<int>& all_positions);
    ThreadingInstructions(const std::vector<std::vector<int>>& starts,
                          const std::vector<std::vector<double>>& tmrcas,
                          const std::vector<std::vector<int>>& targets,
                          const std::vector<std::vector<int>>& mismatches,
                          const std::vector<int>& _positions, int _start, int _end);

public:
    int start = 0;
    int end = 0;
    int num_samples = 0;
    int num_sites = 0;
    std::vector<int> positions;
    std::vector<ThreadingInstruction> instructions;
};

std::tuple<
    std::vector<std::vector<std::vector<int>>>,
    std::vector<std::vector<std::vector<double>>>,
    std::vector<std::vector<std::vector<int>>>,
    std::vector<std::vector<std::vector<int>>>,
    std::vector<std::vector<int>>,
    std::vector<int>,
    std::vector<int>
>
batch_threading_instructions(
    const std::vector<std::vector<int>>& starts,
    const std::vector<std::vector<double>>& tmrcas,
    const std::vector<std::vector<int>>& targets,
    const std::vector<std::vector<int>>& mismatches,
    const std::vector<int>& positions,
    int num_batches
);

#endif // THREADS_ARG_THREADING_INSTRUCTIONS_HPP