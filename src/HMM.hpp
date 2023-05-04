#ifndef DEMOGRAPHY
#define DEMOGRAPHY
#include "Demography.hpp"
#endif DEMOGRAPHY

#include <vector>
#include <iostream>

using std::ostream;

class HMM {
public:
  // HMM data
  int num_states;
  std::vector<double> expected_times;

  std::vector<std::vector<double>> non_transition_score;
  std::vector<std::vector<double>> transition_score;
  std::vector<std::vector<double>> hom_score;
  std::vector<std::vector<double>> het_score;
  void compute_recombination_scores(std::vector<double> cm_sizes);
  void compute_mutation_scores(std::vector<double> bp_sizes, double mutation_rate);
  std::vector<double> compute_expected_times(Demography demography, int K);

  // HMM working quantities
  std::vector<std::vector<double>> trellis;
  std::vector<std::vector<unsigned short>> pointers;

  HMM(Demography demography, std::vector<double> bp_sizes, std::vector<double> cm_sizes, double mutation_rate, int K);
  HMM() {};
  std::vector<int> breakpoints(std::vector<bool> observations, int start);
};
