#ifndef THREADS_ARG_MATCHER_HPP
#define THREADS_ARG_MATCHER_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>

// for a certain interval, store the matches for all samples
class MatchGroup {
public:
  int num_samples;
  std::unordered_map<int, std::unordered_set<int>> match_candidates;
  std::vector<std::unordered_map<int, int>> match_candidates_counts;
  std::vector<std::vector<std::pair<int, int>>> top_four_maps;
  double cm_position;
  MatchGroup(int _num_samples, double cm_position);
  MatchGroup(const std::vector<int>& target_ids,
             const std::vector<std::unordered_set<int>>& matches, const double _cm_position);
  void filter_matches(int min_matches);
  void insert_tops_from(MatchGroup& other);
  void clear();
};

class Matcher {

private:
  int min_matches = 0;
  int sites_processed = 0;
  int next_query_site_idx = 0;
  int match_group_idx = 0;
  int min_match_length = 1;
  std::vector<MatchGroup> match_groups;
  std::vector<int> sorting;
  std::vector<int> next_sorting;
  std::vector<int> permutation;

public:
  int neighborhood_size = 0;
  std::vector<int> query_sites;
  std::vector<int> match_group_sites;
  int num_samples = 0;
  int num_sites = 0;
  double query_interval_size = 0.0;
  // matches in these groups are considered together in the hmm
  double match_group_interval_size = 0.0;
  std::vector<double> genetic_positions;
  Matcher(int _n, const std::vector<double>& _genetic_positions, double _query_interval_size,
          double _match_group_interval_size, int _neighborhood_size, int _min_matches);

  // Do all the work
  void process_site(const std::vector<int>& genotype);
  void propagate_adjacent_matches();
  void clear();
  std::vector<MatchGroup> get_matches();
  std::vector<std::vector<std::unordered_set<int>>>
  serializable_matches(std::vector<int>& target_ids);
  std::vector<double> cm_positions();

  std::vector<int> get_sorting();
  std::vector<int> get_permutation();
};

#endif // THREADS_ARG_MATCHER_HPP
