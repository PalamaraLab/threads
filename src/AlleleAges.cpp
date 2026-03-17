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

#include "AlleleAges.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

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
                // Self-referencing targets don't extend carrier chains
                if (static_cast<size_t>(target) != i) {
                    path_lengths[i] = path_lengths[target] + 1;
                } else {
                    path_lengths[i] = 1;
                }
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
            // Skip self-referencing targets to avoid infinite loops
            if (next_sample == static_cast<int>(start_tmp)) {
                break;
            }
            double this_tmrca = threading_iterators.at(start_tmp).current_tmrca;
            running_max = std::max(running_max, this_tmrca);
            tmrcas.at(next_sample) = running_max;
            start_tmp = next_sample;
        }

        for (int i = 0; i < num_samples; i++) {
            if (tmrcas.at(i) < 0) {
                int target = threading_iterators.at(i).current_target;
                double tmrca = threading_iterators.at(i).current_tmrca;
                if (target == i || target < 0 || tmrcas.at(target) < 0) {
                    // Self-ref or unfilled target: use own tmrca
                    tmrcas.at(i) = tmrca;
                } else {
                    tmrcas.at(i) = std::max(tmrcas.at(target), tmrca);
                }
            }
        }

        // Sort samples by tmrca, then sweep to find the threshold that
        // maximizes: carriers_at_or_below(t) + non_carriers_above(t).
        struct TmrcaSample {
            double tmrca;
            int genotype;
        };
        std::vector<TmrcaSample> sorted_samples;
        sorted_samples.reserve(num_samples);
        int total_non_carriers = 0;
        for (int i = 0; i < num_samples; i++) {
            sorted_samples.push_back({tmrcas[i], genotypes[i]});
            if (genotypes[i] == 0) total_non_carriers++;
        }
        // NaN-safe sort: put NaN values last
        std::sort(sorted_samples.begin(), sorted_samples.end(),
                  [](const TmrcaSample& a, const TmrcaSample& b) {
                      if (std::isnan(a.tmrca)) return false;
                      if (std::isnan(b.tmrca)) return true;
                      return a.tmrca < b.tmrca;
                  });

        // Initial score: threshold below all tmrcas → 0 carriers correct,
        // all non-carriers correct (they're all above threshold).
        int score = total_non_carriers;
        int best_score = score;
        double best_boundary = sorted_samples.front().tmrca;
        double next_boundary = sorted_samples.front().tmrca;

        // Sweep through sorted samples in groups of equal tmrca
        size_t i_sweep = 0;
        size_t n_sorted = sorted_samples.size();
        while (i_sweep < n_sorted) {
            double current_t = sorted_samples[i_sweep].tmrca;
            // Stop at NaN values (they are sorted to the end)
            if (std::isnan(current_t)) break;
            // Process all samples at this tmrca
            int carriers_at_t = 0;
            int non_carriers_at_t = 0;
            size_t group_end = i_sweep;
            while (group_end < n_sorted && sorted_samples[group_end].tmrca == current_t) {
                if (sorted_samples[group_end].genotype > 0)
                    carriers_at_t++;
                else
                    non_carriers_at_t++;
                group_end++;
            }
            // Moving threshold to include this group:
            // carriers become correctly classified (+), non-carriers become incorrect (-)
            score += carriers_at_t - non_carriers_at_t;

            if (score > best_score) {
                best_score = score;
                best_boundary = current_t;
                next_boundary = (group_end < n_sorted)
                    ? sorted_samples[group_end].tmrca
                    : current_t;
            }
            i_sweep = group_end;
        }

        if (best_score > max_score) {
            max_score = best_score;
            allele_age = (best_boundary == next_boundary)
                ? best_boundary + 1
                : (best_boundary + next_boundary) / 2.;
        }
    }

    estimated_ages.push_back(allele_age);
    increment_site();
}

std::vector<double> AgeEstimator::get_inferred_ages() const {
    return estimated_ages;
}
