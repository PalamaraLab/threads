#include "Demography.hpp"
#include <algorithm> // for std::sort
#include <iostream>
#include <math.h>
#include <stdexcept> // to throw errors
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

Demography::Demography(std::vector<double> _sizes, std::vector<double> _times)
    : times(_times), sizes(_sizes), std_times(std::vector<double>()) {
  if (times.size() != sizes.size()) {
    cerr << "Need times and sizes of equal length in demography\n";
    exit(1);
  }
  std::vector<double> deltas;
  for (int i = 0; i < times.size(); i++) {
    if (times[i] < 0 || sizes.size() <= 0) {
      cerr << "Need non-negative times and strictly positive sizes.\n";
      exit(1);
    }
    if (i > 0 && times[i] <= times[i - 1]) {
      cerr << "Demography times must be strictly increasing, found ";
      cerr << times[i] << " after " << times[i - 1] << endl;
      exit(1);
    }
  }

  if (times[0] > 0) {
    cerr << "Demography must start at time 0.0, found " << times[0] << endl;
    exit(1);
  }

  // Compute times in standard coalescent space
  int K = times.size();
  std_times.push_back(0.0);
  for (int i = 1; i < K; i++) {
    double d = (times[i] - times[i - 1]) / sizes[i - 1];
    std_times.push_back(std_times[i - 1] + d);
  }

  // Compute the expected pairwise coalescent time
  // expected_time = expected_branch_length(2);//0;
  expected_time = std_to_gen(1); // 0;
  // for (int i = 0; i < K; i++) {
  //   double T1 = static_cast<double>(times[i]);
  //   double gamma_k = 1. / sizes[i];
  //   if (i < K - 1) {
  //     double T2 = static_cast<double>(times[i + 1]);
  //     expected_time += ((gamma_k * T1 + 1) * std::exp(-gamma_k * T1) - (gamma_k * T2 + 1) *
  //     std::exp(-gamma_k * T2)) / gamma_k;
  //   } else {
  //     expected_time += (gamma_k * T1 + 1) * std::exp(-gamma_k * T1) / gamma_k;
  //   }
  // }
}

double Demography::std_to_gen(const double t) {
  if (t < 0) {
    cerr << "Can only convert non-negative times to std.\n";
    exit(1);
  }
  // Find the highest i s.t. std_times[i] <= t.
  int i =
      std::distance(std_times.begin(), std::upper_bound(std_times.begin(), std_times.end(), t)) - 1;
  return times[i] + (t - std_times[i]) * sizes[i];
}

// /**
//  * @brief Expected length of branch number N
//  *
//  * @param N
//  * @return double
//  */
// double Demography::expected_branch_length(const int N)
// {
//   if (N <= 2) {
//     return expected_time;
//   }
//   int K = times.size();
//   double binom_factor = N * (N - 1) / 2;
//   // Compute the expected pairwise coalescent time
//   double e_time = 0.;
//   for (int i = 0; i < K; i++) {
//     double T1 = static_cast<double>(times[i]);
//     double gamma_k = binom_factor / sizes[i];
//     if (i < K - 1) {
//       double T2 = static_cast<double>(times[i + 1]);
//       e_time += ((gamma_k * T1 + 1) * std::exp(-gamma_k * T1) - (gamma_k * T2 + 1) *
//       std::exp(-gamma_k * T2)) / gamma_k;
//     } else {
//       e_time += (gamma_k * T1 + 1) * std::exp(-gamma_k * T1) / gamma_k;
//     }
//   }
//   return N * e_time;
// }

/**
 * @brief Expected length of branch number N
 *
 * @param N
 * @return double
 */
double Demography::expected_branch_length(const int N) {
  return std_to_gen(2. / N);
  // if (N >= 2) {
  //   return std_to_gen(2 / N);
  // } //else {
  //   return expected_time
  // }
  // ****
  // * rest of this is untested and may have bugs
  // ****
  // int K = times.size();
  // double binom_factor = N * (N + 1) / 2;
  // // Compute the expected pairwise coalescent time
  // double e_time = 0.;
  // for (int k = 0; k < K; k++) {
  //   double T1 = times[k];
  //   double gamma_k = 1 / sizes[k];
  //   double lambda_k = binom_factor * gamma_k;

  //   double coal_fac = std::exp(-binom_factor * std_times[k]);

  //   if (k < K - 1) {
  //     double T2 = times[k + 1];
  //     e_time += coal_fac * ((lambda_k * T1 + 1) - (lambda_k * T2 + 1) * std::exp(-lambda_k * (T2
  //     - T1))) / lambda_k;
  //     // e_time += ((gamma_k * T1 + 1) * std::exp(-gamma_k * T1) - (gamma_k * T2 + 1) *
  //     std::exp(-gamma_k * T2)) / gamma_k;
  //   } else {
  //     e_time += coal_fac * (lambda_k * T1 + 1) / lambda_k;
  //   }
  // }
  // return N * e_time;
}

ostream& operator<<(ostream& os, const Demography& d) {
  for (int i = 0; i < d.sizes.size(); i++) {
    cout << d.times[i] << " " << d.sizes[i] << " " << d.std_times[i] << endl;
  }
  // os << "Node for sample " << node.sample_ID << " carrying allele " << node.genotype;
  return os;
}
