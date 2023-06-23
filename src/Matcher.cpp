#include "Matcher.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

// for a given interval, this contains all the matches for all the samples
MatchGroup::MatchGroup(int _num_samples, double _cm_position)
    : num_samples(_num_samples), cm_position(_cm_position) {
  for (int i = 0; i < num_samples; i++) {
    match_candidates_counts.push_back(std::unordered_map<int, int>()); // could also be a multiset?
  }
}

MatchGroup::MatchGroup(const std::vector<int>& target_ids,
                       const std::vector<std::unordered_set<int>>& matches,
                       const double _cm_position)
    : cm_position(_cm_position) {
  if (target_ids.size() != matches.size()) {
    throw std::runtime_error("Inconsistent target/matches sizes");
  }
  for (int i = 0; i < target_ids.size(); i++) {
    match_candidates[target_ids[i]] = matches[i];
  }
}

void MatchGroup::filter_matches(int min_matches) {
  // First set the candidates for this group
  for (int i = 0; i < num_samples; i++) {
    match_candidates[i] = {};
    if (i < 100) {
      for (int j = 0; j < i; j++) {
        match_candidates.at(i).insert(j);
      }
    }
    else if (i < 1000) {
      for (auto counts : match_candidates_counts.at(i)) {
        if (counts.second >= std::min(2, min_matches)) {
          match_candidates.at(i).insert(counts.first);
        }
      }
    }
    else {
      for (auto counts : match_candidates_counts.at(i)) {
        if (counts.second >= min_matches) {
          match_candidates.at(i).insert(counts.first);
        }
      }
    }
    // special case if nothing was found below threshold (can happen)
    if (match_candidates.at(i).size() == 0) {
      int tmp_min_matches = min_matches;
      while (match_candidates.at(i).size() == 0 && tmp_min_matches > 0) {
        for (auto counts : match_candidates_counts.at(i)) {
          if (counts.second >= tmp_min_matches) {
            match_candidates.at(i).insert(counts.first);
          }
        }
        tmp_min_matches--;
      }
    }
  }

  // Then determine top 4 candidates for neighbouring groups
  top_four_maps.reserve(num_samples);
  for (int i = 0; i < num_samples; i++) {
    top_four_maps.emplace_back(std::min(4, (int) match_candidates.at(i).size()));
    // top_four_maps.resize();
    std::partial_sort_copy(match_candidates_counts.at(i).begin(),
                           match_candidates_counts.at(i).end(), top_four_maps.at(i).begin(),
                           top_four_maps.at(i).end(),
                           [](std::pair<int, int> const& l, std::pair<int, int> const& r) {
                             return l.second > r.second;
                           });
  }
}

void MatchGroup::insert_tops_from(MatchGroup& other) {
  for (int i = 1; i < num_samples; i++) {
    // cout << "tops " << i << endl;
    for (auto p : other.top_four_maps.at(i)) {
      match_candidates.at(i).insert(p.first);
    }
  }
}

Matcher::Matcher(int _n, const std::vector<double>& _genetic_positions, double _query_interval_size,
                 double _match_group_interval_size, int _L, int _min_matches)
    : num_samples(_n), genetic_positions(_genetic_positions),
      query_interval_size(_query_interval_size),
      match_group_interval_size(_match_group_interval_size), L(_L), min_matches(_min_matches) {
  if (genetic_positions.size() <= 2) {
    throw std::runtime_error("Need at least 3 sites, found " +
                             std::to_string(genetic_positions.size()));
  }
  num_sites = genetic_positions.size();
  sites_processed = 0;

  // Check maps are strictly increasing
  for (int i = 0; i < num_sites - 1; i++) {
    if (genetic_positions.at(i + 1) <= genetic_positions.at(i)) {
      std::string prompt = "Genetic coordinates must be strictly increasing, found ";
      prompt += std::to_string(genetic_positions[i + 1]) + " after " +
                std::to_string(genetic_positions[i]);
      throw std::runtime_error(prompt);
    }
  }

  int query_site_idx = 1;
  int match_group_site_idx = 1;
  double gen_pos_offset = genetic_positions[0];
  // for (double cm : genetic_positions) {
  match_group_sites = {0};
  for (int i = 0; i < genetic_positions.size(); i++) {
    double cm = genetic_positions.at(i);
    if (cm > gen_pos_offset + query_interval_size * query_site_idx) {
      query_sites.push_back(i);
      while (cm > gen_pos_offset + query_interval_size * query_site_idx) {
        query_site_idx++;
      }
    }
    if (cm > gen_pos_offset + match_group_interval_size * match_group_site_idx) {
      match_group_sites.push_back(i);
      while (cm > gen_pos_offset + match_group_interval_size * match_group_site_idx) {
        match_group_site_idx++;
      }
    }
  }

  if (query_sites.size() == 0) {
    std::string prompt = "Query interval size larger than genetic map, interval size (cM): ";
    prompt += std::to_string(query_interval_size) + ", map size (cM): ";
    prompt += std::to_string(genetic_positions.back() - genetic_positions.at(0));
    throw std::runtime_error(prompt);
  }
  next_query_site_idx = 0;

  if (match_group_sites.size() == 0) {
    std::string prompt = "No match groups found, interval size (cM): ";
    prompt += std::to_string(match_group_interval_size) + ", map size (cM): ";
    prompt += std::to_string(genetic_positions.back() - genetic_positions.at(0));
    throw std::runtime_error(prompt);
  }
  match_group_idx = 0;
  cout << "found " << query_sites.size() << " query sites and " << match_group_sites.size()
       << " match_group_sites" << endl;

  match_groups.reserve(match_group_sites.size());
  for (int match_group_site : match_group_sites) {
    match_groups.emplace_back(num_samples, genetic_positions[match_group_site]);
  }

  sorting.reserve(num_samples);
  next_sorting.reserve(num_samples);
  permutation.reserve(num_samples);
  for (int i = 0; i < num_samples; i++) {
    sorting.push_back(i);
    next_sorting.push_back(i);
    permutation.push_back(i);
  }
}

void Matcher::process_site(const std::vector<int>& genotype) {
  if (sites_processed >= num_sites) {
    throw std::runtime_error("all sites have already been processed");
  }

  int allele_count = 0;
  for (int g : genotype) {
    if (g == 1) {
      allele_count++;
    }
  }
  int counter0 = 0;
  int counter1 = 0;
  if (genotype.size() != num_samples) {
    throw std::runtime_error("invalid genotype vector size");
  }

  for (int i = 0; i < num_samples; i++) {
    if (genotype.at(sorting.at(i)) == 1) {
      next_sorting[num_samples - allele_count + counter1] = sorting.at(i);
      counter1++;
    }
    else if (genotype.at(sorting.at(i)) == 0) {
      next_sorting[counter0] = sorting.at(i);
      counter0++;
    }
    else {
      std::string prompt = "invalid genotype" + std::to_string(genotype.at(sorting.at(i)));
      throw std::runtime_error(prompt);
    }
  }
  sorting = next_sorting;


  // do the threading-neighbour queries
  if (match_group_idx < match_group_sites.size() - 1 &&
      sites_processed >= match_group_sites.at(match_group_idx + 1)) {
    // this check is awkward, rewrite
    match_group_idx++;
    match_groups.at(match_group_idx - 1).filter_matches(min_matches);
  }

  if (next_query_site_idx < query_sites.size() &&
      sites_processed == query_sites.at(next_query_site_idx)) {
    // do argsort(next_sorting)
    for (int i = 0; i < num_samples; i++) {
      permutation[sorting.at(i)] = i;
    }
    next_query_site_idx++;

    std::set<int> threaded = {permutation.at(0)};

    for (int i = 1; i < num_samples; i++) {
      std::vector<int> matches;
      int allele = genotype.at(i);
      matches.reserve(L);
      auto iter = threaded.insert(permutation.at(i));
      auto iter_up = iter.first;
      auto iter_down = iter.first;
      // check if genotypes are identical, just to be sure
      while (matches.size() < L && (iter_down != threaded.begin() || iter_up != threaded.end())) {
        if (iter_down != threaded.begin()) {
          iter_down--;
          matches.push_back(sorting.at(*iter_down));
        }
        if (matches.size() < L && iter_up != threaded.end()) {
          iter_up++;
          if (iter_up != threaded.end()) {
            matches.push_back(sorting.at(*iter_up));
          }
        }
      }

      for (int m : matches) {
        std::unordered_map<int, int>& mmmap =
            match_groups.at(match_group_idx).match_candidates_counts.at(i);
        if (m >= i) {
          throw std::runtime_error("illegal match candidate " + std::to_string(m) +
                                   ", wth is going on");
        }
        // match_groups.at(match_group_idx).match_candidates.at(i).insert(m);
        if (!mmmap.count(m)) {
          mmmap[m] = 1;
        }
        else {
          mmmap[m]++;
        }
      }
    }

    // special case for last query
    if (next_query_site_idx == query_sites.size()) {
      match_groups.at(match_group_sites.size() - 1).filter_matches(min_matches);
    }
  }
  sites_processed++;
}

// bool Matcher::cached_match(int i, int j) {
//   for (int l = 0; l < cache_size; l++) {
//     if (cached_genotypes.at(l).at(i) != cached_genotypes.at(l).at(j)) {
//       return false;
//     }
//   }
//   return true;
// }

void Matcher::propagate_adjacent_matches() {
  // for (MatchGroup& group : match_groups) {
  //   group.set_candidates(min_matches);
  // }

  for (int i = 1; i < match_groups.size(); i++) {
    MatchGroup& group = match_groups.at(i);
    MatchGroup& prev = match_groups.at(i - 1);
    group.insert_tops_from(prev);
    prev.insert_tops_from(group);
  }
}

std::vector<MatchGroup> Matcher::get_matches() {
  return match_groups;
}

// This returns a list (groups) of lists (targets) of sets (matches)
std::vector<std::vector<std::unordered_set<int>>>
Matcher::serializable_matches(std::vector<int>& target_ids) {
  std::vector<std::vector<std::unordered_set<int>>> serialized_matches;
  serialized_matches.reserve(match_groups.size());
  for (MatchGroup& match_group : match_groups) {
    std::vector<std::unordered_set<int>> current_group_matches;
    for (int target_id : target_ids) {
      current_group_matches.push_back(match_group.match_candidates.at(target_id));
    }
    serialized_matches.push_back(current_group_matches);
  }
  return serialized_matches;
}

std::vector<double> Matcher::cm_positions() {
  std::vector<double> cms;
  cms.reserve(match_groups.size());
  for (MatchGroup& match_group : match_groups) {
    cms.push_back(match_group.cm_position);
  }
  return cms;
}

std::vector<int> Matcher::get_sorting() {
  return sorting;
}

std::vector<int> Matcher::get_permutation() {
  return permutation;
}
