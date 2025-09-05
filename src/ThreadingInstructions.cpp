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

#include "ThreadingInstructions.hpp"
#include "GenotypeIterator.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>


// perf: Note that args are passed by value rather than ref as they run faster
// when used with std::move below.
ThreadingInstruction::ThreadingInstruction(std::vector<int> _starts,
                                           std::vector<double> _tmrcas,
                                           std::vector<int> _targets,
                                           std::vector<int> _mismatches) :
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

ThreadingInstructions::ThreadingInstructions(std::vector<ThreadingInstruction> _instructions, std::vector<int> _positions) :
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

ThreadingInstructions ThreadingInstructions::sub_range(const int range_start, const int range_end) const {
    // Validate range and get corresponding start end index in positions
    if (range_start > range_end) {
        std::ostringstream oss;
        oss << "range_start " << range_start << " cannot be greater than range_end " << range_end;
        throw std::runtime_error(oss.str());
    }

    const auto range_start_it = std::lower_bound(positions.begin(), positions.end(), range_start);
    if ((range_start < 0) || (range_start_it == positions.end())) {
        std::ostringstream oss;
        oss << "range_start " << range_start << " is not a valid position";
        throw std::runtime_error(oss.str());
    }
    const int range_start_idx = range_start_it - positions.begin();

    const auto range_end_it = std::lower_bound(positions.begin(), positions.end(), range_end);
    if ((range_end < 0) || (range_end_it == positions.end())) {
        std::ostringstream oss;
        oss << "range_end " << range_end << " is not a valid position";
        throw std::runtime_error(oss.str());
    }
    const int range_end_idx = range_end_it - positions.begin();

    // Create range positions from indexed subset
    std::vector<int> range_positions{
        positions.begin() + range_start_idx,
        positions.begin() + range_end_idx + 1
    };

    std::vector<ThreadingInstruction> range_instructions;
    for (const auto& instruction : instructions) {
        std::vector<int> sub_starts;
        std::vector<double> sub_tmrcas;
        std::vector<int> sub_targets;
        std::vector<int> sub_mismatches;

        // Create sub-vectors based on range
        const size_t num_starts = instruction.starts.size();
        for (size_t pos_idx = 0; pos_idx < num_starts; ++pos_idx) {
            const int seg_start = instruction.starts[pos_idx];
            const int seg_end = (pos_idx == num_starts - 1) ?
                std::numeric_limits<int>::max() :
                instruction.starts[pos_idx + 1];

            if ((range_start <= seg_end) && (range_end >= seg_start)) {
                sub_starts.push_back(seg_start);
                sub_tmrcas.push_back(instruction.tmrcas[pos_idx]);
                sub_targets.push_back(instruction.targets[pos_idx]);
            }
        }

        // Recompute mismatches based on start indexes within range
        for (const int mismatch_idx : instruction.mismatches) {
            const int pos = positions[mismatch_idx];
            if ((pos >= range_start) && (pos <= range_end)) {
                sub_mismatches.push_back(mismatch_idx - range_start_idx);
            }
        }

        if (!sub_starts.empty() || !sub_mismatches.empty()) {
            range_instructions.emplace_back(
                std::move(sub_starts),
                std::move(sub_tmrcas),
                std::move(sub_targets),
                std::move(sub_mismatches)
            );
        }
    }

    return ThreadingInstructions{
        std::move(range_instructions),
        std::move(range_positions)
    };
}

std::vector<double> ThreadingInstructions::left_multiply(const std::vector<double>& x, bool diploid, bool normalize) {
    // Left-multiplication of the genotype matrix by a vector of doubles

    // Check input vector lengths are correct
    if (diploid) {
        if (x.size() != num_samples / 2) {
            std::ostringstream oss;
            oss << "Input vector must have length " << num_samples / 2 << ".";
            throw std::runtime_error(oss.str());
        }
    } else {
        if (x.size() != num_samples) {
            std::ostringstream oss;
            oss << "Input vector must have length " << num_samples << ".";
            throw std::runtime_error(oss.str());
        }
    }

    // Initialize genotype traversal
    GenotypeIterator gi = GenotypeIterator(*this);
    std::size_t site_counter = 0;
    std::vector<double> out(num_sites);

    while (gi.has_next_genotype()) {
        // Fetch the next genotype
        const std::vector<int>& g = gi.next_genotype();

        // Initialize the next entry
        double entry = 0.0;

        if (normalize) {
            // If we want to normalize, we need the mean and standard deviation of g.
            double ac = 0.0;
            for (auto a : g) {
                ac += a;
            }
            if (diploid) {
                // We do the diploid standard deviation by hand
                double mu = 2.0 * ac / num_samples;
                double sample_var = 0.0;
                for (std::size_t i=0; i < x.size(); i++) {
                    int h = g[2 * i] + g[2 * i + 1];
                    double d = h - mu;
                    sample_var += d * d;
                }
                sample_var /= (num_samples / 2);

                double std = std::sqrt(sample_var);
                for (std::size_t i=0; i < x.size(); i++) {
                    int h = g[2 * i] + g[2 * i + 1];
                    double w = x[i];
                    entry += w * (h - mu) / std;
                }
            } else {
                double mu = ac / num_samples;
                double std = std::sqrt(mu * (1 - mu));
                for (std::size_t i=0; i < g.size(); i++) {
                    double w = x[i];
                    entry += w * (g[i] - mu) / std;
                }
            }
        } else {
            for (std::size_t i=0; i < g.size(); i++) {
                double w = diploid ? x[i / 2] : x[i];
                entry += w * g[i];
            }
        }
        out[site_counter] = entry;
        site_counter++;
    }
    return out;
}

std::vector<double> ThreadingInstructions::right_multiply(const std::vector<double>& x, bool diploid, bool normalize) {
    // Right-multiplication of the genotype matrix by a vector of doubles

    // Check input vector lengths are correct
    if (x.size() != num_sites) {
        std::ostringstream oss;
        oss << "Input vector must have length " << num_samples / 2 << ".";
        throw std::runtime_error(oss.str());
    }

    GenotypeIterator gi = GenotypeIterator(*this);
    std::size_t site_counter = 0;
    if (diploid) {
        // Initialize output
        std::vector<double> out(num_samples / 2, 0.0);
        if (normalize) {
            while (gi.has_next_genotype()) {
                // Fetch the next genotype
                const std::vector<int>& g = gi.next_genotype();

                // If we want to normalize, we need the mean and standard deviation of g.
                double ac = 0.0;
                for (auto a : g) {
                    ac += a;
                }

                // We do the diploid standard deviation by hand
                const double mu = 2.0 * ac / num_samples;
                double sample_var = 0.0;
                for (std::size_t i=0; i < out.size(); i++) {
                    int h = g[2 * i] + g[2 * i + 1];
                    double d = h - mu;
                    sample_var += d * d;
                }
                sample_var /= (num_samples / 2);
                const double std = std::sqrt(sample_var);

                const double w = x[site_counter] / std;
                for (std::size_t i=0; i < out.size(); i++) {
                    const int h = g[2 * i] + g[2 * i + 1];
                    out[i] += w * (h - mu);
                }
                site_counter++;
            }
        } else {
            while (gi.has_next_genotype()) {
                // Fetch the next genotype
                const std::vector<int>& g = gi.next_genotype();
                const double w = x[site_counter];
                for (std::size_t i=0; i < out.size(); i++) {
                    const int h = g[2 * i] + g[2 * i + 1];
                    out[i] += w * h;
                }
                site_counter++;
            }
        }
        return out;
    } else {
        // Initialize output
        std::vector<double> out(num_samples, 0.0);
        if (normalize) {
            while (gi.has_next_genotype()) {
                // Fetch the next genotype
                const std::vector<int>& g = gi.next_genotype();
                double ac = 0.0;
                for (auto a : g) {
                    ac += a;
                }

                // Normalization constants
                double mu = ac / num_samples;
                double std = std::sqrt(mu * (1 - mu));
                const double w = x[site_counter] / std;
                for (std::size_t i=0; i < out.size(); i++) {
                    out[i] += w * (g[i] - mu);
                }
                site_counter++;
            }
        } else {
            while (gi.has_next_genotype()) {
                // Fetch the next genotype
                const std::vector<int>& g = gi.next_genotype();
                const double w = x[site_counter];
                for (std::size_t i=0; i < out.size(); i++) {
                    out[i] += w * g[i];
                }
                site_counter++;
            }
        }
        return out;
    }
}
