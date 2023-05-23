#include <unordered_map>
#include <unordered_set>
#include <vector>

// NB this recycles a lot of code from Matcher.hpp. One day, maybe, try to clean that up.
class ImputationMatcher {

private:
  int sites_processed;
  int next_query_site_idx;

  // pbwt quantities
  std::vector<int> sorting;
  std::vector<int> next_sorting;

  // querying quantities
  std::vector<int> ref_sorting;
  // std::vector<int> divergence;
  // std::vector<int> next_divergence;

public:
  int L;
  std::vector<int> query_sites;
  int num_samples;
  int num_reference;
  int num_target;
  int num_sites;
  double query_interval_size;
  std::unordered_map<int, std::unordered_set<int>> match_sets;
  // neighbourhood size
  int neighborhood_size;
  std::vector<double> genetic_positions;
  ImputationMatcher(int _n_ref, int _n_target, const std::vector<double>& _genetic_positions,
                    double _query_interval_size, int _L);

  // Do all the work
  void process_site(const std::vector<int>& genotype);
  std::unordered_map<int, std::unordered_set<int>> get_matches();
  // std::vector<std::vector<std::unordered_set<int>>>
  // serializable_matches(std::vector<int>& target_ids);

  // std::vector<double> cm_positions();

  std::vector<int> get_sorting();
  // std::vector<int> get_permutation();
  // std::vector<int> get_divergence();
};
