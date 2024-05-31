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

#ifndef THREADS_ARG_DEMOGRAPHY_HPP
#define THREADS_ARG_DEMOGRAPHY_HPP

#include <iostream>
#include <vector>

/// This class is a wrapper for simple coalescence time queries under a piecewise-constant
/// demography
class Demography {
public:
  Demography(std::vector<double> _times, std::vector<double> _sizes);

  /// Map a time in the standard coalescent to generations under this demography
  double std_to_gen(const double t);

  /// The expected branch length of a new branch in a tree with N leaves
  double expected_branch_length(const int N);

  /// Stream output
  friend std::ostream& operator<<(std::ostream& os, const Demography& demography);

public:
  std::vector<double> times;     ///< These are in generations
  std::vector<double> sizes;     ///< These are *haploid*
  std::vector<double> std_times; ///< Normalised coalescence times
  double expected_time = 0.0;    ///< Expected pairwise coalescent time
};

#endif // THREADS_ARG_DEMOGRAPHY_HPP
