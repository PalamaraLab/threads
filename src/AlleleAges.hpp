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

#ifndef THREADS_ARG_ALLELE_AGES_HPP
#define THREADS_ARG_ALLELE_AGES_HPP

#include "ThreadingInstructions.hpp"
#include <iostream>
#include <vector>

class AgeEstimator {
private:
    int sites_processed = 0;
    int num_sites = 0;
    int num_samples = 0;
    std::vector<int> positions;
    std::vector<ThreadingInstructionIterator> threading_iterators;
public:
    AgeEstimator(ThreadingInstructions& instructions);
    void process_site(const std::vector<int>& genotypes);
    void increment_site();
    const std::vector<double> get_inferred_ages();
public:
    std::vector<double> estimated_ages;
};

#endif // THREADS_ARG_ALLELE_AGES_HPP