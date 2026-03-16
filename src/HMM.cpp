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

#include "HMM.hpp"

#include <boost/math/distributions/exponential.hpp>
#include <iostream>
#include <vector>

HMM::HMM(Demography demography, std::vector<double> bp_sizes, std::vector<double> cm_sizes,
         double mutation_rate, int K)
    : num_states(K) {
  expected_times = compute_expected_times(demography, K);

  compute_recombination_scores(cm_sizes);
  compute_mutation_scores(bp_sizes, mutation_rate);

  trellis.reserve(bp_sizes.size());
  pointers.reserve(bp_sizes.size());
  for (std::size_t i = 0; i < bp_sizes.size(); i++) {
    trellis.emplace_back(num_states, 0.0);
    pointers.emplace_back(num_states, 0);
  }
}

std::vector<double> HMM::compute_expected_times(Demography demography, const int K) {
  std::vector<double> result;
  result.reserve(K);
  double k = static_cast<double>(num_states);
  boost::math::exponential e;

  for (int i = 1; i <= K; i++) {
    double t = demography.std_to_gen(quantile(e, (i - 0.5) / k));
    result.push_back(t);
  }
  return result;
}

void HMM::compute_recombination_scores(std::vector<double> cm_sizes) {
  const double log_K = std::log(num_states);
  non_transition_score.reserve(cm_sizes.size());
  transition_score.reserve(cm_sizes.size());
  for (std::size_t i = 0; i < cm_sizes.size(); i++) {
    std::vector<double> non_trans(num_states);
    std::vector<double> trans(num_states);
    for (int k = 0; k < num_states; k++) {
      const double l = 2. * 0.01 * cm_sizes[i] * expected_times[k];
      const double t = std::log1p(-std::exp(-l)) - log_K;

      // log-prob of transitioning
      trans[k] = t;

      // log-prob of *not* transitioning
      non_trans[k] = std::log(std::exp(-l) + std::exp(t));
    }
    transition_score.push_back(std::move(trans));
    non_transition_score.push_back(std::move(non_trans));
  }
}

void HMM::compute_mutation_scores(std::vector<double> bp_sizes, double mutation_rate) {
  hom_score.reserve(bp_sizes.size());
  het_score.reserve(bp_sizes.size());
  for (std::size_t i = 0; i < bp_sizes.size(); i++) {
    std::vector<double> hom(num_states);
    std::vector<double> het(num_states);
    for (int k = 0; k < num_states; k++) {
      // TODO: use mean-bp sizes here as in the main algorithm
      const double l = 2. * mutation_rate * bp_sizes[i] * expected_times[k];

      // log-prob of mutating
      het[k] = std::log1p(-std::exp(-l));

      // log-prob of *not* mutating
      hom[k] = -l;
    }
    het_score.push_back(std::move(het));
    hom_score.push_back(std::move(hom));
  }
}

// PSMC-like algorithm to get breakpoints
// (This is an inefficient quadratic implementation, consider
// longer burn-in and linear-time computation for future versions.)
std::vector<int> HMM::breakpoints(std::vector<bool> observations, int start) {
  // Viterbi
  // Initialize
  int neighborhood_size = static_cast<int>(observations.size());
  std::vector<unsigned short> z(neighborhood_size);
  for (int i = 0; i < num_states; i++) {
    double score = observations[0] ? het_score[start][i] : hom_score[start][i];
    trellis[start][i] = score;
  }

  if (num_states > std::numeric_limits<unsigned short>::max()) {
    throw std::runtime_error("Unable to store breakpoints for more than 2^16 states");
  }

  // Main routine
  double score = 0.0;
  unsigned short running_argmax = 0;
  for (int j = 1; j < neighborhood_size; j++) {
    const int js = j + start;
    for (int i = 0; i < num_states; i++) {
      // Hoist mut_score out of the inner k-loop: it only depends on i, not k
      const double mut_score = observations[j] ? het_score[js][i] : hom_score[js][i];
      double running_max = 0;
      for (int k = 0; k < num_states; k++) {
        double rec_score =
            k == i ? non_transition_score[js][k] : transition_score[js][k];

        score = trellis[j - 1 + start][k] + rec_score + mut_score;

        if (score > running_max || k == 0) {
          running_max = score;
          running_argmax = static_cast<unsigned short>(k);
        }
      }
      trellis[js][i] = running_max;
      pointers[js][i] = running_argmax;
    }
  }

  // Get best path
  double running_max = trellis[start + neighborhood_size - 1][0];
  unsigned short argmax = 0;
  for (int k = 1; k < num_states; k++) {
    double s = trellis[start + neighborhood_size - 1][k];
    if (s > running_max) {
      running_max = s;
      argmax = static_cast<unsigned short>(k);
    }
  }

  // Traceback
  z[neighborhood_size - 1] = argmax;
  for (int j = neighborhood_size - 1; j >= 1; j--) {
    z[j - 1] = pointers[j + start][z[j]];
  }

  // Break it up
  std::vector<int> breakpoints;
  breakpoints.push_back(0 + start);
  for (int j = 1; j < neighborhood_size; j++) {
    if (z[j - 1] != z[j]) {
      breakpoints.push_back(j + start);
    }
  }

  return breakpoints;
}
