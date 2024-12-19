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

#include "ThreadingInstructions.hpp"
#include <limits>
#include <iostream>
#include <vector>


ThreadingInstruction::ThreadingInstruction(const std::vector<int> _starts,
                                           const std::vector<double> _tmrcas,
                                           const std::vector<int> _targets,
                                           const std::vector<int> _mismatches) : 
    starts(_starts), tmrcas(_tmrcas), targets(_targets), mismatches(_mismatches) {
    num_segments = starts.size();
    if (tmrcas.size() != num_segments || targets.size() != num_segments) {
      throw std::runtime_error("Mismatching lengths of threading instruction input");
    }
}

ThreadingInstructionIterator::ThreadingInstructionIterator(ThreadingInstruction& _instruction) :
    instruction(_instruction) {
    if (instruction.num_segments == 0) {
        current_tmrca = 0;
        current_target = -1;
    } else {
        current_tmrca = instruction.tmrcas.at(0);
        current_target = instruction.targets.at(0);
    }
    if (instruction.num_segments <= 1) {
        next_segment_start = std::numeric_limits<int>::max();
    } else {
        next_segment_start = instruction.starts.at(0);
    }
}

void ThreadingInstructionIterator::increment_site(int position) {
    // Search forward and update current target/tmrca for this position
    while (next_segment_start <= position && next_segment_start != std::numeric_limits<int>::max()) {
        current_segment++;
        current_target = instruction.targets.at(current_segment);
        current_tmrca = instruction.tmrcas.at(current_segment);
        if (current_segment < instruction.num_segments - 1) {
            next_segment_start = instruction.starts.at(current_segment + 1);
        } else {
            next_segment_start = std::numeric_limits<int>::max();
        }
    }
    sites_processed++;
}

ThreadingInstructions::ThreadingInstructions(const std::vector<ThreadingInstruction>& _instructions, const std::vector<int>& _positions) :
    instructions(_instructions), positions(_positions) {
    start = positions.at(0);
    end = positions.back() + 1;
    num_samples = instructions.size();
    num_sites = positions.size();
}

ThreadingInstructions::ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int _start, const int _end, const std::vector<int>& all_positions) :
    start(_start), end(_end) {
    num_samples = paths.size();

    for (auto pos : all_positions) {
        if (start <= pos && pos < end) {
            positions.push_back(pos);
        } 
        if (pos > end) {
            break;
        }
    }
    num_sites = positions.size();

    std::vector<int> starts;
    std::vector<int> targets;
    std::vector<double> tmrcas;
    for (auto path : paths) {
        std::vector<int> mismatches;
        path.map_positions(all_positions);
        std::tie(starts, targets, tmrcas) = path.dump_data_in_range(start, end);
        for (auto het_idx : path.het_sites) {
            int het_pos = all_positions.at(het_idx);
            if ((start <= het_pos) && (het_pos < end)) {
                mismatches.push_back(het_idx);
            }
        }
        instructions.push_back(ThreadingInstruction(starts, tmrcas, targets, mismatches));
    }
}

ThreadingInstructions::ThreadingInstructions(const std::vector<std::vector<int>>& starts,
                        const std::vector<std::vector<double>>& tmrcas,
                        const std::vector<std::vector<int>>& targets,
                        const std::vector<std::vector<int>>& mismatches,
                        const std::vector<int>& _positions, int _start, int _end) :
    positions(_positions), start(_start), end(_end) {
    num_samples = starts.size();
    num_sites = positions.size();
    for (int i = 0; i < num_samples; i++) {
        instructions.push_back(ThreadingInstruction(starts.at(i), tmrcas.at(i), targets.at(i), mismatches.at(i)));
    }
}

std::vector<std::vector<int>> ThreadingInstructions::all_starts() {
    std::vector<std::vector<int>> out;
    for (auto& instruction : instructions) {
        out.push_back(instruction.starts);
    }
    return out;
}

std::vector<std::vector<double>> ThreadingInstructions::all_tmrcas() {
    std::vector<std::vector<double>> out;
    for (auto& instruction : instructions) {
        out.push_back(instruction.tmrcas);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_targets() {
    std::vector<std::vector<int>> out;
    for (auto& instruction : instructions) {
        out.push_back(instruction.targets);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_mismatches() {
    std::vector<std::vector<int>> out;
    for (auto& instruction : instructions) {
        out.push_back(instruction.mismatches);
    }
    return out;
}