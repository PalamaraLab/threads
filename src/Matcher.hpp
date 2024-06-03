// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef THREADS_ARG_MATCHER_HPP
#define THREADS_ARG_MATCHER_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>

/// for a certain interval, store the matches for all samples
class MatchGroup {
public:
  MatchGroup(int _num_samples, double cm_position);
  MatchGroup(const std::vector<int>& target_ids,
             const std::vector<std::unordered_set<int>>& matches, const double _cm_position);
  void filter_matches(int min_matches);
  void insert_tops_from(MatchGroup& other);
  void clear();

public:
  int num_samples = 0;
  std::unordered_map<int, std::unordered_set<int>> match_candidates;
  std::vector<std::unordered_map<int, int>> match_candidates_counts;
  std::vector<std::vector<std::pair<int, int>>> top_four_maps;
  double cm_position = 0.0;
};

class Matcher {
public:
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

public:
  int num_samples = 0;
  std::vector<double> genetic_positions;
  double query_interval_size = 0.0;
  double match_group_interval_size = 0.0;
  int neighborhood_size = 0;
  std::vector<int> query_sites;
  std::vector<int> match_group_sites;
  int num_sites = 0;
  // matches in these groups are considered together in the hmm

private:
  int min_matches = 0;
  int sites_processed = 0;
  int next_query_site_idx = 0;
  int match_group_idx = 0;
  std::vector<MatchGroup> match_groups;
  std::vector<int> sorting;
  std::vector<int> next_sorting;
  std::vector<int> permutation;
};

#endif // THREADS_ARG_MATCHER_HPP
