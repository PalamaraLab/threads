#include "Demography.hpp"
#include <vector>
#include <iostream>
#include <algorithm> // for std::sort
#include <stdexcept> // to throw errors

using std::endl;
using std::cout;
using std::cerr;

Demography::Demography(std::vector<double> _sizes, std::vector<double> _times) :
  times(_times), sizes(_sizes), std_times(std::vector<double>()) 
{
  if (times.size() != sizes.size()) {
    cerr << "Need times and sizes of equal length in demography\n";
    exit(1);
  }
  for (int i = 0; i < times.size(); i++) {
    if (times[i] < 0 || sizes.size() <= 0) {
      cerr << "Need non-negative times and strictly positive sizes.\n";
      exit(1);
    }
    if (i > 0 && times[i] <= times[i - 1]) {
      cerr << "Times must be strictly increasing.\n";
      exit(1);
    }
  }

  if (times[0] > 0) {
    cerr << "Demography must start at time 0.0, found " << times[0] << endl;
    exit(1);
  }

  std_times.push_back(0.0);
  for (int i = 1; i < times.size(); i++) {
    double d = (times[i] - times[i - 1]) / sizes[i - 1];
    std_times.push_back(std_times[i - 1] + d);
  }
}

double Demography::std_to_gen(double t) 
{
  if (t < 0) {
    cerr << "Can only convert non-negative times to std.\n";
    exit(1);
  }
  // Find the highest i s.t. std_times[i] <= t.
  int i = std::distance(std_times.begin(), std::upper_bound(std_times.begin(), std_times.end(), t)) - 1;
  return times[i] + (t - std_times[i]) * sizes[i];
}

ostream& operator<<(ostream& os, const Demography& d) {
  for (int i = 0; i < d.sizes.size(); i++) {
    cout << d.times[i] << " " << d.sizes[i] << " " << d.std_times[i] << endl;
  }
  // os << "Node for sample " << node.sample_ID << " carrying allele " << node.genotype;
  return os;
}
