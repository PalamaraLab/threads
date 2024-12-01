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

#include "DataConsistency.hpp"
#include <limits>
#include <iostream>
#include <vector>


InstructionConverter::InstructionConverter(ThreadingInstruction _instructions, std::size_t _instruction_index, std::size_t start_position) :
    instructions(_instructions), instruction_index(_instruction_index) {
    num_segments = instructions.num_segments;
    converted_segment_start = start_position;
    if (instruction_index == 0) {
        current_target = -1;
    } else {
        current_target = instructions.targets.at(0);
    }
    sites_processed = 0;
    if (instructions.num_segments > 1) {
        next_segment_start = instructions.starts.at(1);
    } else {
        next_segment_start = std::numeric_limits<int>::max();
    }
}

void InstructionConverter::break_segment(double new_lower_bound, double new_upper_bound, int position, int new_target) {
    double threads_tmrca = instructions.tmrcas.at(current_segment);

    if (threads_tmrca >= current_upper_bound) {
        // Just below the current upper bound
        double epsilon = std::min(current_upper_bound / 1000, (current_upper_bound - current_lower_bound) / 2);
        threads_tmrca = current_upper_bound - epsilon;
    } else if (threads_tmrca <= current_lower_bound) {
        // Just above the current lower bound
        double epsilon = std::min(current_lower_bound / 1000, (current_upper_bound - current_lower_bound) / 2);
        threads_tmrca = current_lower_bound + epsilon;
    }

    if ((threads_tmrca <= current_lower_bound) || (current_upper_bound <= threads_tmrca)) {
        throw std::runtime_error("Threading time exceeds allowed bounds, this is not allowed.");
    }

    new_tmrcas.push_back(threads_tmrca);
    new_starts.push_back(converted_segment_start);
    new_targets.push_back(current_target);
    converted_segment_start = position;
    current_target = new_target;
    current_upper_bound = new_upper_bound;
    current_lower_bound = new_lower_bound;
}

void InstructionConverter::increment_site(std::size_t position) {
    // Check if we need to update current_bounds
    if (position == next_segment_start) {
        // Entering new threading segment, we need to wrap up current segment and start a new one
        break_segment(0, std::numeric_limits<double>::max(), position, current_target);
        current_segment++;
        current_target = instructions.targets.at(current_segment);
        if (current_segment < num_segments - 1) {
            next_segment_start = instructions.starts.at(current_segment + 1);
        } else {
            next_segment_start = std::numeric_limits<int>::max();
        }
    }
}

void InstructionConverter::update_target(double new_lower_bound, double new_upper_bound, std::size_t position, int new_target) {
    if (new_target != current_target) {
        if (position != converted_segment_start) {
            break_segment(new_lower_bound, new_upper_bound, position, new_target);
        } else {
            current_target = new_target;
        }
    }
}

void InstructionConverter::evaluate_bounds(std::vector<int>& genotypes, std::size_t position, double allele_age) {
    if (instruction_index == 0) {
        sites_processed++;
        return;
    }

    // extract relevant genotypes
    int gt_focal = genotypes.at(instruction_index);
    int gt_target = genotypes.at(current_target);

    // Get allele age bounds for this site, if any
    double local_upper_bound = std::numeric_limits<double>::max();
    double local_lower_bound = 0;

    // If 1/1, both are carriers and must coalsce below the mutation.
    // If 0/1 or 1/0, only one is a carrier and these lineages must coalesce
    // above the mutation.
    // If 0/0, this mutation does not provide any information on this coalescence event.
    if (gt_focal == 1 && gt_target == 1) {
        local_upper_bound = allele_age;
    } else if (gt_focal == 1 || gt_target == 1) {
        local_lower_bound = allele_age;
    }

    if ((local_upper_bound <= current_lower_bound) || (local_lower_bound >= current_upper_bound)) {
        // Bounds don't match, we need to wrap up current segment and start a new one
        break_segment(local_lower_bound, local_upper_bound, position, current_target);
    } else {
        // Otherwise, we update current bounds and continue
        current_lower_bound = std::max(local_lower_bound, current_lower_bound);
        current_upper_bound = std::min(local_upper_bound, current_upper_bound);
    }
    sites_processed++;
}

ThreadingInstruction InstructionConverter::parse_converted_instructions() {
    if (instruction_index == 0) {
        return instructions;
    }
    if (current_segment != instructions.num_segments - 1) {
        std::cout << current_segment << " should be " << instructions.num_segments - 1 << "\n";
        throw std::runtime_error("Haven't processed all segments in a set of instructions");
    }
    // Wrapping up, so we break up and save the current segment, setting
    // The "-1"s are unused placeholders.
    break_segment(-1, -1, -1, -1);
    return ThreadingInstruction(new_starts, new_tmrcas, new_targets, new_mismatches);
}

ConsistencyWrapper::ConsistencyWrapper(const std::vector<std::vector<int>>& starts,
                                       const std::vector<std::vector<double>>& tmrcas,
                                       const std::vector<std::vector<int>>& targets,
                                       const std::vector<std::vector<int>>& mismatches,
                                       const std::vector<int>& _physical_positions,
                                       const std::vector<double>& _allele_ages)
    : allele_ages(_allele_ages), physical_positions(_physical_positions) {
    num_sites = allele_ages.size();
    num_samples = starts.size();
    if (num_sites == 0) {
        throw std::runtime_error("Found 0 sites");
    }
    if (num_samples == 0) {
        throw std::runtime_error("Found 0 samples");
    }
    if (tmrcas.size() != num_samples || targets.size() != num_samples || mismatches.size() != num_samples) {
      throw std::runtime_error("Mismatching lengths of threading instruction input");
    }
    if (physical_positions.size() != num_sites) {
        throw std::runtime_error("Sites don't match ages");
    }

    for (std::size_t i = 0; i < num_samples; i++) {
        ThreadingInstruction instructions = ThreadingInstruction(starts.at(i), tmrcas.at(i), targets.at(i), mismatches.at(i));
        InstructionConverter converter = InstructionConverter(instructions, i, physical_positions.at(0));
        instruction_converters.push_back(converter);
    }

    std::cout << "Will convert " << num_samples << " threading instructions across " << num_sites << " sites.\n";
}

void ConsistencyWrapper::process_site(std::vector<int>& genotypes) {
    double allele_age = allele_ages.at(sites_processed);
    // An age of 0 shouldn't happen, but just in case, we interpret it as just really small
    if (allele_age <= 0) {
        allele_age = 1e-3;
    }
    std::size_t position = physical_positions.at(sites_processed);
    int counter = 0;
    int first_carrier = -1;
    for (InstructionConverter& converter : instruction_converters) {
        converter.increment_site(position);
        int new_target = converter.current_target;
        if (genotypes.at(counter) == 1) {
            // The threading instructions of carriers must be a carrier-only graph rooted
            // at the first carrier
            if (first_carrier == -1) {
                first_carrier = counter;
            } else {
                // Otherwise we try traversing the local threading graph to find another carrier
                int current_target = converter.current_target;
                while (current_target != -1 && genotypes.at(current_target) != 1) {
                    current_target = instruction_converters.at(current_target).current_target;
                }
                if (current_target > 0 && genotypes.at(current_target) == 1) {
                    new_target = current_target;
                } else {
                    // If that doesn't work just use the "first" carrier
                    new_target = first_carrier;
                }
            }
        }
        if (new_target != converter.current_target) {
            // Force a new threading segment bounded by [0, allele_age) and
            // using the new target
            converter.update_target(0., allele_age, position, new_target);
        }
        converter.evaluate_bounds(genotypes, position, allele_age);

        if (counter > 0 && genotypes.at(converter.current_target) != genotypes.at(counter)) {
            converter.new_mismatches.push_back(counter);
        }
        counter++;
    }
    sites_processed++;
}

ThreadingInstructions ConsistencyWrapper::get_consistent_instructions() {
    if (sites_processed != num_sites) {
        throw std::runtime_error("Haven't processed all sites yet, cannot output converted instructions");
    }

    // Make output threading instructions
    std::vector<ThreadingInstruction> output_instructions;
    for (InstructionConverter converter : instruction_converters) {
        output_instructions.push_back(converter.parse_converted_instructions());
    }

    return ThreadingInstructions(output_instructions, physical_positions);
}
