#include <vector>
#include <iostream>
using std::ostream;

class Demography {
private:
  std::vector<double> std_times;

public:
  // These are in generations
  std::vector<double> times;
  // These are *haploid*
  std::vector<double> sizes;

  Demography(std::vector<double> _times, std::vector<double> _sizes);

  // Map a time in the standard coalescent to generations under this demography
  double std_to_gen(double t);
  
  // Output
  friend ostream& operator<<(ostream& os, const Demography& demography);
};
