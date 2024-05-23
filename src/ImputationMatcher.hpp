#ifndef THREADS_INFER_IMPUTATION_MATCHER_HPP
#define THREADS_INFER_IMPUTATION_MATCHER_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>

// NB this recycles a lot of code from Matcher.hpp
class ImputationMatcher {

private:
  int sites_processed;
  int next_query_site_idx;

  // pbwt quantities
  std::vector<int> sorting;
  std::vector<int> next_sorting;

  // querying quantities
  std::vector<int> ref_sorting;

public:
  // This is a version of the Threads haplotype matching algorithm
  // adapted to be used for imputation.
  // TODO: include a second pass through data here to get divergence values and not do that using
  // Threads-fastLS
  int L;
  std::vector<int> query_sites;
  int num_samples;
  int num_reference;
  int num_target;
  int num_sites;
  double query_interval_size;
  std::unordered_map<int, std::unordered_set<int>> match_sets;
  int neighborhood_size;
  std::vector<double> genetic_positions;
  ImputationMatcher(int _n_ref, int _n_target, const std::vector<double>& _genetic_positions,
                    double _query_interval_size, int _L);

  // Do all the work
  void process_site(const std::vector<int>& genotype);
  std::unordered_map<int, std::unordered_set<int>> get_matches();

  std::vector<int> get_sorting();
};

#endif // THREADS_INFER_IMPUTATION_MATCHER_HPP
