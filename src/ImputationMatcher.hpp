#ifndef THREADS_ARG_IMPUTATION_MATCHER_HPP
#define THREADS_ARG_IMPUTATION_MATCHER_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>

/// This is a version of the Threads haplotype matching algorithm
/// adapted to be used for imputation.
/// NB this recycles a lot of code from Matcher.hpp
class ImputationMatcher {
public:
  ImputationMatcher(int _n_ref, int _n_target, const std::vector<double>& _genetic_positions,
                    double _query_interval_size, int _neighborhood_size);

  void process_site(const std::vector<int>& genotype);
  const std::unordered_map<int, std::unordered_set<int>>& get_matches() const;
  const std::vector<int>& get_sorting() const;

public:
  // TODO: include a second pass through data here to get divergence values and not do that using
  // Threads-fastLS
  int neighborhood_size = 0;
  std::vector<int> query_sites;
  int num_samples = 0;
  int num_reference = 0;
  int num_target = 0;
  int num_sites = 0;
  double query_interval_size = 0.0;
  std::unordered_map<int, std::unordered_set<int>> match_sets;
  std::vector<double> genetic_positions;

private:
  int sites_processed = 0;
  int next_query_site_idx = 0;

  // pbwt quantities
  std::vector<int> sorting;
  std::vector<int> next_sorting;

  // querying quantities
  std::vector<int> ref_sorting;
};

#endif // THREADS_ARG_IMPUTATION_MATCHER_HPP
