// This file is part of the Threads software suite.
// Copyright (C) 2025 Threads Developers.
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

#include "GenotypeIterator.hpp"

#include <vector>

GenotypeIterator::GenotypeIterator(const ThreadingInstructions& instructions) {
    positions = instructions.positions;
    for (auto& instruction : instructions.instructions) {
        threading_iterators.emplace_back(instruction, positions);
    }
    num_samples = instructions.num_samples;
    out_genotype = std::vector(num_samples, 0);
    current_site = 0;
    num_sites = instructions.num_sites;
    current_position = positions.front();
    reference_genome = std::vector(instructions.num_sites, 0);
    for (auto i : instructions.instructions.front().mismatches) {
        reference_genome[i] = 1;
    }
}

const std::vector<int>& GenotypeIterator::next_genotype() {
    // Recover the current site
    for (int i = 0; i < num_samples; i++) {
        if (i == 0) {
            out_genotype[i] = reference_genome.at(current_site);
        } else {
            auto& iterator = threading_iterators.at(i);
            out_genotype[i] = iterator.is_mismatch ? 1 - out_genotype.at(iterator.current_target) : out_genotype.at(iterator.current_target);
        }
    }

    // Increment the iterators
    current_site++;
    if (current_site < positions.size()) {
        for (auto& iterator : threading_iterators) {
            iterator.increment_site(positions.at(current_site));
        }
        current_position = positions.at(current_site);
    } else {
        current_position = -1;
    }

    return out_genotype;
}

bool GenotypeIterator::has_next_genotype() {
    return current_site < num_sites;
}
