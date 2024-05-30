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

  // FIXME Alex review with Arni, worth reserving up front or using emplace_back?
  for (std::size_t i = 0; i < bp_sizes.size(); i++) {
    std::vector<double> trellis_row(num_states, 0.0);
    std::vector<unsigned short> pointer_row(num_states, 0);
    trellis.push_back(trellis_row);
    pointers.push_back(pointer_row);
  }
}

std::vector<double> HMM::compute_expected_times(Demography demography, const int K) {
  std::vector<double> result;
  double k = static_cast<double>(num_states);
  boost::math::exponential e;

  for (int i = 1; i <= K; i++) {
    double t = demography.std_to_gen(quantile(e, (i - 0.5) / k));
    result.push_back(t);
  }
  return result;
}

void HMM::compute_recombination_scores(std::vector<double> cm_sizes) {
  for (std::size_t i = 0; i < cm_sizes.size(); i++) {
    non_transition_score.push_back(std::vector<double>());
    transition_score.push_back(std::vector<double>());
    for (int k = 0; k < num_states; k++) {
      double t = expected_times[k];
      const double l = 2. * 0.01 * cm_sizes[i] * t;
      const double trans = std::log1p(-std::exp(-l)) - std::log(num_states);

      // log-prob of transitioning
      transition_score[i].push_back(trans);

      // log-prob of *not* transitioning
      non_transition_score[i].push_back(std::log(std::exp(-l) + std::exp(trans)));
    }
  }
}

void HMM::compute_mutation_scores(std::vector<double> bp_sizes, double mutation_rate) {
  for (std::size_t i = 0; i < bp_sizes.size(); i++) {
    hom_score.push_back(std::vector<double>());
    het_score.push_back(std::vector<double>());
    for (int k = 0; k < num_states; k++) {
      double t = expected_times[k];

      // TODO: use mean-bp sizes here as in the main algorithm
      const double l = 2. * mutation_rate * bp_sizes[i] * t;

      // log-prob of mutating
      het_score[i].push_back(std::log1p(-std::exp(-l)));

      // log-prob of *not* mutating
      hom_score[i].push_back(-l);
    }
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

  // FIXME review with Arni, defensive check for assignment to 'k' below
  if (num_states > std::numeric_limits<unsigned short>::max()) {
    throw std::runtime_error("Unable to store breakpoints for more than 2^16 states");
  }

  // Main routine
  double score = 0.0;
  unsigned short running_argmax = 0;
  for (int j = 1; j < neighborhood_size; j++) {
    for (int i = 0; i < num_states; i++) {
      double running_max = 0;
      for (int k = 0; k < num_states; k++) {
        double mut_score = observations[j] ? het_score[j + start][i] : hom_score[j + start][i];
        double rec_score =
            k == i ? non_transition_score[j + start][k] : transition_score[j + start][k];

        score = trellis[j - 1 + start][k] + rec_score + mut_score;

        if (score > running_max || k == 0) {
          running_max = score;
          running_argmax = static_cast<unsigned short>(k);
        }
      }
      trellis[j + start][i] = running_max;
      pointers[j + start][i] = running_argmax;
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
