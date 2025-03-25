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

#ifndef THREADS_ARG_GENOTYPE_RECOVERY_HPP
#define THREADS_ARG_GENOTYPE_RECOVERY_HPP

#include "ThreadingInstructions.hpp"

#include <iostream>
#include <vector>

class GenotypeIterator {
// This is a wrapper to loop over and reconstruct genotypes from Threading instructions
private:
    std::vector<ThreadingInstructionIterator> threading_iterators;
    std::vector<int> out_genotype;
    std::vector<int> reference_genome;

public:
    // Constructor
    GenotypeIterator(ThreadingInstructions& _instructions);

    // Simple iterator functions
    std::vector<int>& next_genotype();
    bool has_next_genotype();

public:
    std::vector<int> positions;
    int num_sites = 0;
    int num_samples = 0;
    int current_site = 0;
    int current_position = 0;
};

#endif // THREADS_ARG_GENOTYPE_RECOVERY_HPP
