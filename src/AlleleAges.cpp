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

#include "AlleleAges.hpp"

#include <numeric>
#include <execution>
#include <vector>
#include <unordered_map>
#include <map>

#include <boost/container/flat_set.hpp>

AgeEstimator::AgeEstimator(const ThreadingInstructions& instructions) {
    num_samples = instructions.num_samples;
    positions = instructions.positions;
    threading_iterators.reserve(instructions.instructions.size());
    for (auto& instruction : instructions.instructions) {
        threading_iterators.emplace_back(instruction, positions);
    }
}

void AgeEstimator::increment_site() {
    int position = positions.at(sites_processed);
    for (auto& iterator : threading_iterators) {
        iterator.increment_site(position);
    }
    sites_processed++;
}

void AgeEstimator::process_site(const std::vector<int>& genotypes) {
    // Find carrier clusters
    // Algorithm: find paths in the marginal threading tree consisting only of carriers
    // Index those by "start" (highest index)
    std::vector<size_t> max_path_starts;
    size_t max_path_len = 0;

    std::vector<size_t> path_lengths(genotypes.size(), 0);
    for (size_t i = 0; i < genotypes.size(); i++) {
        if (i == 0) {
            path_lengths[i] = genotypes.at(i);
        } else {
            if (genotypes.at(i) == 1) {
                int target = threading_iterators.at(i).current_target;
                path_lengths[i] = path_lengths[target] + 1;
            } else {
                path_lengths[i] = 0;
            }
        }

        // Find the largest carrier cluster(s)
        // Extract the longest such path(s)
        if (path_lengths[i] == max_path_len) {
            max_path_starts.push_back(i);
        } else if (path_lengths[i] > max_path_len) {
            max_path_len = path_lengths[i];
            max_path_starts = {i};
        }
    }

    int max_score = 0;
    double allele_age = 0;
    for (size_t path_start : max_path_starts) {
        // For each longest path, trace the start node to the root and compute
        // tmrca with each other sample

        // This map keeps track of all tmrcas with the sample at "path_start"
        std::vector<double> tmrcas(num_samples, -1);
        tmrcas.at(path_start) = 0;
        double running_max = 0;
        size_t start_tmp = path_start;
        while (start_tmp > 0) {
            int next_sample = threading_iterators.at(start_tmp).current_target;
            double this_tmrca = threading_iterators.at(start_tmp).current_tmrca;
            running_max = std::max(running_max, this_tmrca);
            tmrcas.at(next_sample) = running_max;
            start_tmp = next_sample;
        }
        if (start_tmp != 0) {
            throw std::runtime_error("Invalid threading instruction traversal.");
        }

        for (int i = 0; i < num_samples; i++) {
            if (tmrcas.at(i) < 0) {
                int target = threading_iterators.at(i).current_target;
                double tmrca = threading_iterators.at(i).current_tmrca;
                tmrcas.at(i) = std::max(tmrcas.at(target), tmrca);
            }
        }

        // Create a sorted unique list of coalescence times as transform_reduce
        // below must be done in order. For performance, boost's flat_set is
        // faster than std::set or sorting a std::vector in this instance.
        boost::container::flat_set<double> unique_tmrcas(tmrcas.begin(), tmrcas.end());

        // For each sample, check its tmrca with path_start and
        // update the score for each tmrca bin accordingly
        std::map<double, int> scores;
        for (double t : unique_tmrcas) {
            scores[t] = std::transform_reduce(
                tmrcas.begin(),
                tmrcas.end(),
                genotypes.begin(),
                0,
                std::plus<>(),
                [t](double tmrca, int genotype) {
                    bool is_carrier = genotype > 0;
                    if (is_carrier && tmrca <= t) {
                        return 1;
                    }
                    if (!is_carrier && tmrca > t) {
                        return 1;
                    }
                    return 0;
                }
            );
        }

        std::vector<double> age_bin_boundaries;
        for (auto const& imap: scores)
            age_bin_boundaries.push_back(imap.first);
        std::vector<double> age_bins;

        for (size_t k = 0; k < age_bin_boundaries.size(); k++) {
            int score = scores.at(age_bin_boundaries.at(k));
            if (score > max_score) {
                max_score = score;
                allele_age = (k == age_bin_boundaries.size() - 1)
                    ? age_bin_boundaries.at(k) + 1
                    : (age_bin_boundaries.at(k) + age_bin_boundaries.at(k + 1)) / 2.;
            }
        }
    }

    estimated_ages.push_back(allele_age);
    increment_site();
}

std::vector<double> AgeEstimator::get_inferred_ages() const {
    return estimated_ages;
}
