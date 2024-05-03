#include "ImputationMatcher.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


ImputationMatcher::ImputationMatcher(int _n_ref, int _n_target,
                                     const std::vector<double>& _genetic_positions,
                                     double _query_interval_size, int _L)
    : num_reference(_n_ref), num_target(_n_target), genetic_positions(_genetic_positions),
      query_interval_size(_query_interval_size), L(_L) {
  if (genetic_positions.size() <= 2) {
    throw std::runtime_error("Need at least 3 sites, found " +
                             std::to_string(genetic_positions.size()));
  }
  num_sites = genetic_positions.size();
  num_samples = num_reference + num_target;
  sites_processed = 0;

  // Check maps are strictly increasing
  // for (int i = 0; i < num_sites - 1; i++) {
  //   if (genetic_positions.at(i + 1) <= genetic_positions.at(i)) {
  //     std::string prompt = "Genetic coordinates must be strictly increasing, found ";
  //     prompt += std::to_string(genetic_positions[i + 1]) + " after " +
  //               std::to_string(genetic_positions[i]);
  //     throw std::runtime_error(prompt);
  //   }
  // }

  int query_site_idx = 1;
  double gen_pos_offset = genetic_positions[0];
  for (int i = 0; i < genetic_positions.size(); i++) {
    double cm = genetic_positions.at(i);
    if (cm > gen_pos_offset + query_interval_size * query_site_idx) {
      query_sites.push_back(i);
      while (cm > gen_pos_offset + query_interval_size * query_site_idx) {
        query_site_idx++;
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

  sorting.reserve(num_samples);
  next_sorting.reserve(num_samples);
  for (int i = 0; i < num_samples; i++) {
    sorting.push_back(i);
    next_sorting.push_back(i);
  }

  ref_sorting.reserve(num_reference);
  for (int i = 0; i < num_reference; i++) {
    ref_sorting.push_back(i);
  }

  match_sets.reserve(num_target);
  for (int i = num_reference; i < num_samples; i++) {
    match_sets[i] = std::unordered_set<int>();
  }
}

void ImputationMatcher::process_site(const std::vector<int>& genotype) {
  if (sites_processed >= num_sites) {
    throw std::runtime_error("all sites have already been processed");
  }

  int allele_count = 0;
  for (int g : genotype) {
    if (g == 1) {
      allele_count++;
    }
  }

  if (genotype.size() != num_samples) {
    throw std::runtime_error("invalid genotype vector size");
  }

  if (next_query_site_idx < query_sites.size() &&
      sites_processed == query_sites.at(next_query_site_idx)) {
    next_query_site_idx++;
    int ref_allele_count = 0;
    for (int i = 0; i < num_reference; i++) {
      if (genotype.at(i) == 1) {
        ref_allele_count++;
      }
    }
    // do matching here
    int ref_counter0 = 0;
    int ref_counter1 = 0;
    std::unordered_map<int, int> target_sort;

    for (int i = 0; i < num_samples; i++) {
      int sample_id = sorting.at(i);
      if (sample_id < num_reference) {
        if (genotype.at(sample_id) == 1) {
          ref_sorting[num_reference - ref_allele_count + ref_counter1] = sample_id;
          ref_counter1++;
        }
        else if (genotype.at(sample_id) == 0) {
          ref_sorting[ref_counter0] = sample_id;
          ref_counter0++;
        }
        else {
          std::string prompt = "invalid genotype" + std::to_string(genotype.at(sample_id));
          throw std::runtime_error(prompt);
        }
      }
      else {
        // Here we do abnormal target-position-finding
        if (genotype.at(sample_id) == 1) {
          target_sort[sample_id] = num_reference - ref_allele_count + ref_counter1;
        }
        else if (genotype.at(sample_id) == 0) {
          target_sort[sample_id] = ref_counter0;
        }
        else {
          std::string prompt = "invalid genotype" + std::to_string(genotype.at(sample_id));
          throw std::runtime_error(prompt);
        }
      }
    }

    for (auto& sorting_twople : target_sort) {
      // get L-sized neighbourhood per target sample around target_sort[target_id] in ref_sorting
      int target_id = sorting_twople.first;
      int sorting_idx = sorting_twople.second;
      int insert_start;
      int insert_end;
      if (sorting_idx < 4) {
        insert_start = 0;
        insert_end = 8;
      }
      else if (sorting_idx > num_reference - 4) {
        insert_start = num_reference - 8;
        insert_end = num_reference;
      }
      else {
        insert_start = sorting_idx - 4;
        insert_end = sorting_idx + 4;
      }
      for (int k = insert_start; k < insert_end; k++) {
        match_sets[target_id].insert(ref_sorting.at(k));
      }
    }
  }

  int counter0 = 0;
  int counter1 = 0;
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
  sites_processed++;
}

std::unordered_map<int, std::unordered_set<int>> ImputationMatcher::get_matches() {
  return match_sets;
}

std::vector<int> ImputationMatcher::get_sorting() {
  return sorting;
}
