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

#ifndef THREADS_ARG_DATA_CONSISTENCY_HPP
#define THREADS_ARG_DATA_CONSISTENCY_HPP

#include "ThreadingInstructions.hpp"
#include <iostream>
#include <limits>
#include <vector>

class InstructionConverter {
private:
    std::vector<int> new_starts;
    std::vector<double> new_tmrcas;
    std::vector<int> new_targets;
    void break_segment(double new_lower_bound, double new_upper_bound, int position, int new_target);
public:
    std::vector<int> new_mismatches;
    InstructionConverter(ThreadingInstruction instructions, std::size_t _instruction_index, std::size_t start_position);
    void process_site(std::vector<int>& genotypes, std::size_t position, double allele_age);
    void increment_site(std::size_t position);
    void update_target(double new_lower_bound, double new_upper_bound, std::size_t position, int new_target);
    void evaluate_bounds(std::vector<int>& genotypes, std::size_t position, double allele_age);

    ThreadingInstruction parse_converted_instructions();
public:
    std::size_t sites_processed = 0;
    std::size_t num_segments = 0; // Number of segments in this set of instructions
    ThreadingInstruction instructions; // The instructions object
    std::size_t instruction_index = 0;  // Index of these instructions in the big set of all instructions
    std::size_t current_segment = 0;  // Index of the current segment in the data traversal
    std::size_t converted_segment_start = 0;
    std::size_t next_segment_start = 0;
    double current_lower_bound = 0;
    double current_upper_bound = std::numeric_limits<double>::max();
    int current_target = 0;
};

class ConsistencyWrapper {
public:
    void process_site(std::vector<int>& genotypes);
    ThreadingInstructions parse_converted_instructions();
    ConsistencyWrapper(const std::vector<std::vector<int>>& starts,
                       const std::vector<std::vector<double>>& tmrcas,
                       const std::vector<std::vector<int>>& targets,
                       const std::vector<std::vector<int>>& mismatches,
                       const std::vector<int>& physical_positions,
                       const std::vector<double>& allele_ages);
    ConsistencyWrapper(ThreadingInstructions& instructions, const std::vector<double>& allele_ages) :
        ConsistencyWrapper(instructions.all_starts(), 
                           instructions.all_tmrcas(),
                           instructions.all_targets(),
                           instructions.all_mismatches(),
                           instructions.positions,
                           allele_ages) {}
    ThreadingInstructions get_consistent_instructions();

public:
    std::vector<double> allele_ages;
    std::vector<int> physical_positions;
    std::size_t sites_processed = 0;
    std::size_t num_sites = 0;
    std::size_t num_samples = 0;
    std::vector<InstructionConverter> instruction_converters;
};

#endif // THREADS_ARG_DATA_CONSISTENCY_HPP