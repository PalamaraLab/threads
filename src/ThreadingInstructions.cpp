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

#include <iostream>
#include <limits>
#include <vector>


ThreadingInstruction::ThreadingInstruction(const std::vector<int>& _starts,
                                           const std::vector<double>& _tmrcas,
                                           const std::vector<int>& _targets,
                                           const std::vector<int>& _mismatches) :
    starts(_starts), tmrcas(_tmrcas), targets(_targets), mismatches(_mismatches) {
    num_segments = starts.size();
    num_mismatches = mismatches.size();
    if (tmrcas.size() != num_segments || targets.size() != num_segments) {
      throw std::runtime_error("Mismatching lengths of threading instruction input");
    }
}

ThreadingInstructionIterator::ThreadingInstructionIterator(
    const ThreadingInstruction& _instruction,
    const std::vector<int>& _positions
) :
    instruction(_instruction), positions(_positions) {
    if (instruction.num_segments == 0) {
        current_tmrca = 0;
        current_target = -1;
    } else {
        current_tmrca = instruction.tmrcas.front();
        current_target = instruction.targets.front();
    }
    if (instruction.num_segments <= 1) {
        next_segment_start = std::numeric_limits<int>::max();
    } else {
        next_segment_start = instruction.starts.at(1);
    }

    if (instruction.num_mismatches == 0) {
        next_mismatch = std::numeric_limits<int>::max();
        is_mismatch = false;
    } else {
        next_mismatch = positions.at(instruction.mismatches.front());
        is_mismatch = instruction.mismatches.front() == 0;
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

    // Search forward and update current mismatch status
    while (next_mismatch < position ) {
        current_mismatch++;
        if (current_mismatch <= instruction.num_mismatches - 1) {
            next_mismatch = positions.at(instruction.mismatches.at(current_mismatch));
        } else {
            next_mismatch = std::numeric_limits<int>::max();
        }
    }
    is_mismatch = position == next_mismatch;

    sites_processed++;
}

ThreadingInstructions::ThreadingInstructions(const std::vector<ThreadingInstruction>& _instructions, const std::vector<int>& _positions) :
    instructions(_instructions), positions(_positions) {
    start = positions.front();
    end = positions.back() + 1;
    num_samples = instructions.size();
    num_sites = positions.size();
}

ThreadingInstructions::ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int _start, const int _end, const std::vector<int>& all_positions) :
    start(_start), end(_end) {
    num_samples = paths.size();

    int pos_idx_offset = -1;
    int i = 0;
    for (auto pos : all_positions) {
        if (start <= pos && pos < end) {
            positions.push_back(pos);
            pos_idx_offset = pos_idx_offset < 0 ? i : pos_idx_offset;
        }
        if (pos > end) {
            break;
        }
        i++;
    }
    pos_idx_offset = pos_idx_offset < 0 ? 0 : pos_idx_offset;
    num_sites = positions.size();

    std::vector<int> starts;
    std::vector<int> targets;
    std::vector<double> tmrcas;
    instructions.reserve(num_samples);
    for (auto path : paths) {
        std::vector<int> mismatches;
        path.map_positions(all_positions);
        std::tie(starts, targets, tmrcas) = path.dump_data_in_range(start, end);
        for (auto het_idx : path.het_sites) {
            int het_pos = all_positions.at(het_idx);
            if ((start <= het_pos) && (het_pos < end)) {
                mismatches.push_back(het_idx - pos_idx_offset);
            }
        }
        instructions.emplace_back(starts, tmrcas, targets, mismatches);
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
    instructions.reserve(num_samples);
    for (int i = 0; i < num_samples; i++) {
        instructions.emplace_back(starts.at(i), tmrcas.at(i), targets.at(i), mismatches.at(i));
    }
}

std::vector<std::vector<int>> ThreadingInstructions::all_starts() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.starts);
    }
    return out;
}

std::vector<std::vector<double>> ThreadingInstructions::all_tmrcas() {
    std::vector<std::vector<double>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.tmrcas);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_targets() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.targets);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_mismatches() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.mismatches);
    }
    return out;
}

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
    const int num_batches
) {
    std::vector<std::vector<std::vector<int>>> batched_starts;
    std::vector<std::vector<std::vector<double>>> batched_tmrcas;
    std::vector<std::vector<std::vector<int>>> batched_targets;
    std::vector<std::vector<std::vector<int>>> batched_mismatches;
    std::vector<std::vector<int>> batched_positions;
    std::vector<int> batched_position_starts;
    std::vector<int> batched_position_ends;

    if (tmrcas.size() != starts.size()) {
        std::ostringstream oss;
        oss << "tmrcas size " << tmrcas.size();
        oss << " does not match starts size " << starts.size();
        throw std::runtime_error(oss.str());
    }

    if (targets.size() != starts.size()) {
        std::ostringstream oss;
        oss << "targets size " << targets.size();
        oss << " does not match starts size " << starts.size();
        throw std::runtime_error(oss.str());
    }

    if (num_batches <= 1) {
        std::ostringstream oss;
        oss << "num_batches must be greater than 1, got " << num_batches;
        throw std::runtime_error(oss.str());
    }

    const int split_size = static_cast<int>(positions.size()) / num_batches;
    int range_start_idx = 0;

    // Step over the positions in segments based on split size
    while (range_start_idx < positions.size()) {
        // Trim any error in position end to actual size
        const int range_end_idx = std::min(
            range_start_idx + split_size - 1,
            static_cast<int>(positions.size())
        );
        if (range_start_idx >= range_end_idx) {
            break;
        }

        // Make positions for this batch as subset of those in range
        const int range_start = positions[range_start_idx];
        const int range_end = positions[range_end_idx];
        std::vector<int> range_positions{
            positions.begin() + range_start_idx,
            positions.begin() + range_end_idx + 1
        };

        std::vector<std::vector<int>> range_starts;
        std::vector<std::vector<double>> range_tmrcas;
        std::vector<std::vector<int>> range_targets;
        std::vector<std::vector<int>> range_mismatches;

        const size_t num_insts = starts.size();
        for (size_t inst_idx = 0; inst_idx < num_insts; ++inst_idx) {
            std::vector<int> sub_starts;
            std::vector<double> sub_tmrcas;
            std::vector<int> sub_targets;
            std::vector<int> sub_mismatches;

            if (tmrcas[inst_idx].size() != starts[inst_idx].size()) {
                std::ostringstream oss;
                oss << "tmrcas size " << tmrcas[inst_idx].size();
                oss << " does not match starts size " << starts[inst_idx].size();
                oss << " for list entry " << inst_idx;
                throw std::runtime_error(oss.str());
            }

            if (targets[inst_idx].size() != starts[inst_idx].size()) {
                std::ostringstream oss;
                oss << "targets size " << targets[inst_idx].size();
                oss << " does not match starts size " << starts[inst_idx].size();
                oss << " for list entry " << inst_idx;
                throw std::runtime_error(oss.str());
            }

            // Create sub-vectors based on range
            const size_t num_starts = starts[inst_idx].size();
            for (size_t pos_idx = 0; pos_idx < num_starts; ++pos_idx) {
                const int seg_start = starts[inst_idx][pos_idx];
                const int seg_end = (pos_idx == num_starts - 1) ?
                    std::numeric_limits<int>::max() :
                    starts[inst_idx][pos_idx + 1];

                if ((range_start <= seg_end) && (range_end >= seg_start)) {
                    sub_starts.push_back(seg_start);
                    sub_tmrcas.push_back(tmrcas[inst_idx][pos_idx]);
                    sub_targets.push_back(targets[inst_idx][pos_idx]);
                }
            }

            // Recompute mismatches based on start indexes within range
            for (const int mismatch_idx : mismatches[inst_idx]) {
                const int pos = positions[mismatch_idx];
                if ((pos >= range_start) && (pos <= range_end)) {
                    sub_mismatches.push_back(mismatch_idx - range_start_idx);
                }
            }

            range_starts.push_back(std::move(sub_starts));
            range_tmrcas.push_back(std::move(sub_tmrcas));
            range_targets.push_back(std::move(sub_targets));
            range_mismatches.push_back(std::move(sub_mismatches));
        }

        batched_starts.push_back(std::move(range_starts));
        batched_tmrcas.push_back(std::move(range_tmrcas));
        batched_targets.push_back(std::move(range_targets));
        batched_mismatches.push_back(std::move(range_mismatches));
        batched_positions.push_back(std::move(range_positions));
        batched_position_starts.push_back(std::move(range_start));
        batched_position_ends.push_back(std::move(range_end));

        range_start_idx += split_size;
    }

    return std::make_tuple(
        batched_starts,
        batched_tmrcas,
        batched_targets,
        batched_mismatches,
        batched_positions,
        batched_position_starts,
        batched_position_ends
    );
}
