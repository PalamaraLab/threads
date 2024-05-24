#ifndef THREADS_ARG_DEMOGRAPHY_HPP
#define THREADS_ARG_DEMOGRAPHY_HPP

#include <vector>
#include <iostream>

// This class is a wrapper for simple coalescence time queries under a piecewise-constant demography
class Demography {
public:
  // These are in generations
  std::vector<double> times;
  // These are *haploid*
  std::vector<double> sizes;
  // Normalised coalescence times
  std::vector<double> std_times;
  // Expected pairwise coalescent time
  double expected_time = 0.0;

  Demography(std::vector<double> _times, std::vector<double> _sizes);

  // Map a time in the standard coalescent to generations under this demography
  double std_to_gen(const double t);
  // The expected branch length of a new branch in a tree with N leaves
  double expected_branch_length(const int N);

  // Output
  friend std::ostream& operator<<(std::ostream& os, const Demography& demography);
};

#endif // THREADS_ARG_DEMOGRAPHY_HPP
