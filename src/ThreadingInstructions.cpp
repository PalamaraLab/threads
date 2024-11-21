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
    current_segment = 0;
    sites_processed = 0;
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

ThreadingInstructions::ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int _start, const int _end, const std::vector<int>& _positions) :
    start(_start), end(_end) {
    num_samples = paths.size();

    for (auto pos : _positions) {
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
        path.map_positions(positions);
        std::tie(starts, targets, tmrcas) = path.dump_data_in_range(start, end);
        for (auto het_idx : path.het_sites) {
            int het_pos = positions.at(het_idx);
            if ((start <= het_pos) && (het_pos < end)) { 
                mismatches.push_back(het_idx);
            }
        }
        instructions.push_back(ThreadingInstruction(starts, tmrcas, targets, mismatches));
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