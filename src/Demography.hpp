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
