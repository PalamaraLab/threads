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
#include <iostream>
#include <vector>


class ThreadingInstruction {
// This is a container class for the instructions inferred for a single sequence
public:
    std::vector<int> starts;
    std::vector<double> tmrcas;
    std::vector<int> targets;
    std::vector<int> mismatches;
    std::size_t num_segments;
    ThreadingInstruction(
        const std::vector<int> _starts,
        const std::vector<double> _tmrcas,
        const std::vector<int> _targets,
        const std::vector<int> _mismatches);
};

class ThreadingInstructionIterator {
// This is a container for iterating through a ThreadingInstruction object site-by-site
private:
    ThreadingInstruction& instruction;
    int sites_processed;
    int next_segment_start;
public:
    int num_sites;
    int current_segment;
    int current_target;
    double current_tmrca;
public:
    void increment_site(const int pos);
    ThreadingInstructionIterator(ThreadingInstruction& _instruction);
};

class ThreadingInstructions {
// This is a container class for threading instructions inferred by Threads
public:
    std::vector<std::vector<int>> all_starts();
    std::vector<std::vector<double>> all_tmrcas();
    std::vector<std::vector<int>> all_targets();
    std::vector<std::vector<int>> all_mismatches();
    std::vector<int> positions;
    ThreadingInstructions(const std::vector<ThreadingInstruction>& _instructions, const std::vector<int>& _positions);
    ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int start, const int end, const std::vector<int>& _positions);
public:
    int start;
    int end;
    int num_samples;
    int num_sites;
    std::vector<ThreadingInstruction> instructions;
};
#endif // THREADS_ARG_THREADING_INSTRUCTIONS_HPP
