#include <vector>
#include <iostream>
using std::ostream;

class Demography {
public:
  // These are in generations
  std::vector<double> times;
  // These are *haploid*
  std::vector<double> sizes;
  // Normalised coalescence times
  std::vector<double> std_times;
  // Expected pairwise coalescent time
  double expected_time;

  Demography(std::vector<double> _times, std::vector<double> _sizes);

  // Map a time in the standard coalescent to generations under this demography
  double std_to_gen(const double t);
  // The expected branch length of a new branch in a tree with N leaves
  double expected_branch_length(const int N);
  
  // Output
  friend ostream& operator<<(ostream& os, const Demography& demography);
};
