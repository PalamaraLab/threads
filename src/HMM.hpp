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
