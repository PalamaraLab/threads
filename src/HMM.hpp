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

#include <iostream>
#include <vector>

/// This class contains the PSMC-like algorithm used to break up segments for small-N inference
class HMM {
public:
  HMM(Demography demography, std::vector<double> bp_sizes, std::vector<double> cm_sizes,
      double mutation_rate, int K);
  HMM() = default;

  void compute_recombination_scores(std::vector<double> cm_sizes);
  void compute_mutation_scores(std::vector<double> bp_sizes, double mutation_rate);
  std::vector<double> compute_expected_times(Demography demography, int K);
  std::vector<int> breakpoints(std::vector<bool> observations, int start);

public:
  // HMM data
  int num_states = 0;
  std::vector<double> expected_times;

  // HMM working quantities
  std::vector<std::vector<double>> trellis;
  std::vector<std::vector<unsigned short>> pointers;

  std::vector<std::vector<double>> non_transition_score;
  std::vector<std::vector<double>> transition_score;
  std::vector<std::vector<double>> hom_score;
  std::vector<std::vector<double>> het_score;
};
