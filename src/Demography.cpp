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

#include "Demography.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <vector>

Demography::Demography(std::vector<double> _sizes, std::vector<double> _times)
    : times(_times), sizes(_sizes), std_times(std::vector<double>()) {
  if (times.size() != sizes.size()) {
    throw std::runtime_error("Demography times and sizes must have equal length");
  }

  for (std::size_t i = 0; i < times.size(); i++) {
    if ((times[i] < 0.0) || (sizes[i] <= 0.0)) {
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
  std::size_t K = times.size();
  std_times.push_back(0.0);
  for (std::size_t i = 1; i < K; i++) {
    double d = (times[i] - times[i - 1]) / sizes[i - 1];
    std_times.push_back(std_times[i - 1] + d);
  }

  // Compute the expected pairwise coalescent time
  expected_time = std_to_gen(1);
}

double Demography::std_to_gen(const double t) {
  if (t < std_times.front()) {
    throw std::runtime_error("Demography can only convert times greater than the first entry");
  }

  // Find the highest i s.t. std_times[i] <= t.
  const auto it = std::upper_bound(std_times.begin(), std_times.end(), t);

  // Defensive check as the t < std_times.front check above should mean `it` is never first
  if (it == std_times.begin()) {
    throw std::runtime_error("Unexpected std_to_gen upper bound finding first element");
  }

  std::size_t i = static_cast<std::size_t>(it - std_times.begin() - 1);
  return times[i] + (t - std_times[i]) * sizes[i];
}

double Demography::expected_branch_length(const int N) {
  return std_to_gen(2. / N);
}

std::ostream& operator<<(std::ostream& os, const Demography& d) {
  for (std::size_t i = 0; i < d.sizes.size(); i++) {
    std::cout << d.times[i] << " " << d.sizes[i] << " " << d.std_times[i] << std::endl;
  }
  return os;
}
