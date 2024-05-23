#include "Demography.hpp"
#include <algorithm> // for std::sort
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept> // to throw errors
#include <vector>


Demography::Demography(std::vector<double> _sizes, std::vector<double> _times)
    : times(_times), sizes(_sizes), std_times(std::vector<double>()) {
  if (times.size() != sizes.size()) {
    throw std::runtime_error("Demography times and sizes must have equal length");
  }

  for (int i = 0; i < times.size(); i++) {
    if (times[i] < 0 || sizes[i] <= 0) {
      throw std::runtime_error("Demography expects non-negative times and strictly positive sizes");
    }

    if (i > 0 && times[i] <= times[i - 1]) {
      std::ostringstream oss;
      oss << "Demography times must be strictly increasing. Found ";
      oss << times[i] << " after " << times[i - 1] << " at index " << i;
      throw std::runtime_error(oss.str());
    }
  }

  if (times[0] > 0) {
    std::ostringstream oss;
    oss << "Demography must start at time 0.0, found " << times[0];
    throw std::runtime_error(oss.str());
  }

  // Compute times in standard coalescent space
  int K = times.size();
  std_times.push_back(0.0);
  for (int i = 1; i < K; i++) {
    double d = (times[i] - times[i - 1]) / sizes[i - 1];
    std_times.push_back(std_times[i - 1] + d);
  }

  // Compute the expected pairwise coalescent time
  expected_time = std_to_gen(1);
}

double Demography::std_to_gen(const double t) {
  if (t < 0) {
    throw std::runtime_error("Demography can only convert non-negative times to std");
  }
  // Find the highest i s.t. std_times[i] <= t.
  int i =
      std::distance(std_times.begin(), std::upper_bound(std_times.begin(), std_times.end(), t)) - 1;
  return times[i] + (t - std_times[i]) * sizes[i];
}

/**
 * @brief Compute the expected length of the N-th branch
 *
 * @param N
 * @return double
 */
double Demography::expected_branch_length(const int N) {
  return std_to_gen(2. / N);
}

std::ostream& operator<<(std::ostream& os, const Demography& d) {
  for (int i = 0; i < d.sizes.size(); i++) {
    std::cout << d.times[i] << " " << d.sizes[i] << " " << d.std_times[i] << std::endl;
  }
  return os;
}
