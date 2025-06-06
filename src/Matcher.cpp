// This file is part of the Threads software suite.
// Copyright (C) 2024-2025 Threads Developers.
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

#include "Matcher.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// For a given interval, this contains all the matches for all the samples
MatchGroup::MatchGroup(int _num_samples, double _cm_position)
    : num_samples(_num_samples), cm_position(_cm_position) {
  for (int i = 0; i < num_samples; i++) {
    match_candidates_counts.push_back(std::unordered_map<int, int>());
  }
}

MatchGroup::MatchGroup(const std::vector<int>& target_ids,
                       const std::vector<std::unordered_set<int>>& matches,
                       const double _cm_position)
    : cm_position(_cm_position) {
  if (target_ids.size() != matches.size()) {
    throw std::runtime_error("Inconsistent target/matches sizes");
  }
  for (int i = 0; i < static_cast<int>(target_ids.size()); i++) {
    match_candidates[target_ids[i]] = matches[i];
  }
}

void MatchGroup::filter_matches(int min_matches) {
  // First set the candidates for this group
  // For sequences of high index, we lower the number of sequences used to save time and memory
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
    else if (i < 10000) {
      for (auto counts : match_candidates_counts.at(i)) {
        if (counts.second >= min_matches) {
          match_candidates.at(i).insert(counts.first);
        }
      }
    }
    else {
      // Don't want too much stuff for very big studies
      for (auto counts : match_candidates_counts.at(i)) {
        if (counts.second >= 2 * min_matches) {
          match_candidates.at(i).insert(counts.first);
        }
      }
    }
    // Special case if nothing was found below threshold (unlikely, but can happen)
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

  // Then determine top 4 candidates for neighboring groups
  top_four_maps.reserve(num_samples);
  for (int i = 0; i < num_samples; i++) {
    top_four_maps.emplace_back(std::min(4, static_cast<int>(match_candidates.at(i).size())));
    std::partial_sort_copy(match_candidates_counts.at(i).begin(),
                           match_candidates_counts.at(i).end(), top_four_maps.at(i).begin(),
                           top_four_maps.at(i).end(),
                           [](std::pair<int, int> const& l, std::pair<int, int> const& r) {
                             return l.second > r.second;
                           });
    match_candidates_counts.at(i).clear();
  }
}

void MatchGroup::insert_tops_from(MatchGroup& other) {
  for (int i = 1; i < num_samples; i++) {
    for (auto p : other.top_four_maps.at(i)) {
      match_candidates.at(i).insert(p.first);
    }
  }
}

void MatchGroup::clear() {
  match_candidates.clear();
  match_candidates_counts.clear();
  top_four_maps.clear();
}

Matcher::Matcher(int _n, const std::vector<double>& _genetic_positions, double _query_interval_size,
                 double _match_group_interval_size, int _neighborhood_size, int _min_matches)
    : num_samples(_n), genetic_positions(_genetic_positions),
      query_interval_size(_query_interval_size),
      match_group_interval_size(_match_group_interval_size), neighborhood_size(_neighborhood_size),
      min_matches(_min_matches) {
  if (genetic_positions.size() <= 2) {
    throw std::runtime_error("Need at least 3 sites, found " +
                             std::to_string(genetic_positions.size()));
  }
  num_sites = static_cast<int>(genetic_positions.size());
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

  // Initialize query sites and match group sites
  int query_site_idx = 1;
  int match_group_site_idx = 1;
  double gen_pos_offset = genetic_positions[0];
  match_group_sites = {0};
  for (int i = 0; i < static_cast<int>(genetic_positions.size()); i++) {
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
  std::cout << "Will use " << query_sites.size() << " query sites and " << match_group_sites.size()
            << " match_group_sites" << std::endl;

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
  // Pass genotypes for a single site through the matcher

  if (sites_processed >= num_sites) {
    throw std::runtime_error("all sites have already been processed");
  }

  // Get allele count
  int allele_count = 0;
  for (int g : genotype) {
    if (g == 1) {
      allele_count++;
    }
  }
  int counter0 = 0;
  int counter1 = 0;
  if (static_cast<int>(genotype.size()) != num_samples) {
    throw std::runtime_error("invalid genotype vector size");
  }

  // PBWT step
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

  // Threading-neighbor queries
  if (match_group_idx < (static_cast<int>(match_group_sites.size()) - 1) &&
      (sites_processed >= match_group_sites.at(match_group_idx + 1))) {
    match_group_idx++;
    match_groups.at(match_group_idx - 1).filter_matches(min_matches);
  }

  // If we've reached a query site, query
  if (next_query_site_idx < static_cast<int>(query_sites.size()) &&
      sites_processed == query_sites.at(next_query_site_idx)) {
    // Get the arg-sort of the sorting
    for (int i = 0; i < num_samples; i++) {
      permutation[sorting.at(i)] = i;
    }
    next_query_site_idx++;

    // Initialize the red-black tree
    std::set<int> threaded = {permutation.at(0)};

    // Insert sequences and query in order
    for (int i = 1; i < num_samples; i++) {
      std::vector<int> matches;
      matches.reserve(neighborhood_size);
      auto iter = threaded.insert(permutation.at(i));
      auto iter_up = iter.first;
      auto iter_down = iter.first;
      // Check if genotypes are identical, just to be sure
      while ((static_cast<int>(matches.size()) < neighborhood_size) &&
             (iter_down != threaded.begin() || iter_up != threaded.end())) {
        if (iter_down != threaded.begin()) {
          iter_down--;
          matches.push_back(sorting.at(*iter_down));
        }
        if (static_cast<int>(matches.size()) < neighborhood_size && iter_up != threaded.end()) {
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
          throw std::runtime_error("Illegal match candidate " + std::to_string(m) +
                                   ", something is very wrong");
        }
        if (!mmmap.count(m)) {
          mmmap[m] = 1;
        }
        else {
          mmmap[m]++;
        }
      }
    }

    // Special case for last query
    if (next_query_site_idx == static_cast<int>(query_sites.size())) {
      match_groups.at(match_group_sites.size() - 1).filter_matches(min_matches);
    }
  }
  sites_processed++;
}

// Propagate top 4 matches from left and right match groups
void Matcher::propagate_adjacent_matches() {
  for (int i = 1; i < static_cast<int>(match_groups.size()); i++) {
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
  std::vector<std::vector<std::unordered_set<int>>> serialized_matches(match_groups.size());
  int group_counter = 0;
  for (MatchGroup& match_group : match_groups) {
    std::vector<std::unordered_set<int>> current_group_matches(target_ids.size());
    int match_counter = 0;
    for (int target_id : target_ids) {
      current_group_matches[match_counter] = std::move(match_group.match_candidates.at(target_id));
      match_group.match_candidates.at(target_id).clear();
      match_counter++;
    }
    serialized_matches[group_counter] = std::move(current_group_matches);
    group_counter++;
  }
  return serialized_matches;
}

void Matcher::clear() {
  match_groups.clear();
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
