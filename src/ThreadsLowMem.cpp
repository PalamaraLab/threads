#include "ThreadsLowMem.hpp"

#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

ThreadsLowMem::ThreadsLowMem(const std::vector<int> _target_ids,
                             const std::vector<double>& _physical_positions,
                             const std::vector<double>& _genetic_positions, std::vector<double> ne,
                             std::vector<double> ne_times, double _mutation_rate, bool _sparse)
    : target_ids(_target_ids), mutation_rate(_mutation_rate),
      physical_positions(_physical_positions), genetic_positions(_genetic_positions),
      sparse(_sparse), demography(Demography(ne, ne_times)) {
  num_samples = static_cast<int>(target_ids.size());
  if (physical_positions.size() != genetic_positions.size()) {
    throw std::runtime_error("Map lengths don't match.");
  }
  else if (physical_positions.size() <= 2) {
    throw std::runtime_error("Need at least 3 sites, found " +
                             std::to_string(physical_positions.size()));
  }
  num_sites = static_cast<int>(physical_positions.size());

  // Check maps are strictly increasing
  for (int i = 0; i < num_sites - 1; i++) {
    if (physical_positions.at(i + 1) < physical_positions.at(i)) {
      std::string prompt = "Physical positions must be increasing, found ";
      prompt += std::to_string(physical_positions[i + 1]) + " after " +
                std::to_string(physical_positions[i]);
      throw std::runtime_error(prompt);
    }
    else if (physical_positions.at(i + 1) == physical_positions.at(i)) {
      physical_positions.at(i + 1) = physical_positions.at(i) + 1;
    }
    if (genetic_positions.at(i + 1) < genetic_positions.at(i)) {
      std::string prompt = "Genetic coordinates must be increasing, found ";
      prompt += std::to_string(genetic_positions[i + 1]) + " after " +
                std::to_string(genetic_positions[i]);
      throw std::runtime_error(prompt);
    }
    else if (genetic_positions.at(i + 1) == genetic_positions.at(i)) {
      genetic_positions.at(i + 1) = genetic_positions.at(i) + 0.0000001;
    }
  }

  // Mean interval size in base-pairs
  mean_bp_size =
      (physical_positions.back() - physical_positions[0]) / static_cast<double>(num_sites - 1);
  for (int target_id : target_ids) {
    segment_indices[target_id] = 0;
    expected_branch_lengths[target_id] = demography.expected_branch_length(target_id + 1);
  }

  // Site counters
  hmm_sites_processed = 0;
  het_sites_processed = 0;
  std::tie(bp_boundaries, bp_sizes) = ThreadsFastLS::site_sizes(physical_positions);

  for (std::size_t i = 0; i < genetic_positions.size(); i++) {
    if (i == genetic_positions.size() - 1) {
      // TO allow for closely spaced markers
      cm_sizes.push_back(0.0000001);
    }
    else {
      cm_sizes.push_back(genetic_positions.at(i + 1) - genetic_positions.at(i));
    }
  }

  // Initialize the psmc-like segment-breaking algorithm
  int min_target_id = *(std::min_element(target_ids.begin(), target_ids.end()));
  if (min_target_id < n_hmm_samples) {
    psmc = HMM(demography, bp_sizes, cm_sizes, mutation_rate, 64);
  }
  else {
    psmc = HMM();
  }
}

// Initialize Threads-Viterbi instances using IDs and genetic positions for each match group
void ThreadsLowMem::initialize_viterbi(std::vector<std::vector<std::unordered_set<int>>>& match_ids,
                                       const std::vector<double>& cm_positions) {
  if (match_ids.size() != cm_positions.size() || match_ids.size() < 1) {
    throw std::runtime_error("Match-data is missing or does not have same shape as genetic map");
  }
  match_groups.reserve(match_ids.size());
  for (std::size_t i = 0; i < match_ids.size(); i++) {
    match_groups.emplace_back(target_ids, match_ids.at(i), cm_positions.at(i));
  }

  match_group_idx = 0;
  hmm_sites_processed = 0;
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    std::vector<int> sample_ids(match_groups.at(0).match_candidates.at(target_id).begin(),
                                match_groups.at(0).match_candidates.at(target_id).end());
    hmms.emplace(target_id, ViterbiState(target_id, sample_ids));
  }
}

// Pass genotypes for a single site through the intialized Threads-Viterbi instances
void ThreadsLowMem::process_site_viterbi(const std::vector<int>& genotype) {
  bool group_change = false;

  if (match_group_idx < (static_cast<int>(match_groups.size()) - 1) &&
      (genetic_positions.at(hmm_sites_processed) >=
       match_groups.at(match_group_idx + 1).cm_position)) {
    match_group_idx++;
    group_change = true;
  }
  double k = 2. * 0.01 * cm_sizes.at(hmm_sites_processed);
  double l = 2. * mutation_rate * bp_sizes.at(hmm_sites_processed);
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    if (group_change) {
      hmms.at(target_id).set_samples(
          match_groups.at(match_group_idx).match_candidates.at(target_id));
    }

    double t = expected_branch_lengths.at(target_id);
    double rho_c = k * t;
    double rho = sparse ? -std::log1p(-std::exp(-(k * t)))
                        : -(std::log1p(-std::exp(-(k * t))) - std::log(target_id));
    double mu_c = l * t;
    double mu = -std::log1p(-std::exp(-(l * t)));
    hmms.at(target_id).process_site(genotype, rho, rho_c, mu, mu_c);
  }
  hmm_sites_processed++;
  return;
}

void ThreadsLowMem::traceback() {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      paths.emplace(target_id, ViterbiPath(0));
    }
    else {

      paths.emplace(target_id, hmms.at(target_id).traceback());
    }
  }
  hmms.clear();
}

void ThreadsLowMem::process_site_hets(const std::vector<int>& genotype) {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      if (genotype.at(0) == 1) {
        paths.at(0).het_sites.push_back(het_sites_processed);
      }
    }
    else {
      ViterbiPath& path = paths.at(target_id);
      int current_seg_idx = segment_indices.at(target_id);
      while (current_seg_idx < (static_cast<int>(path.segment_starts.size()) - 1) &&
             (het_sites_processed >= path.segment_starts.at(current_seg_idx + 1))) {
        current_seg_idx++;
      }
      segment_indices.at(target_id) = current_seg_idx;
      int sample = path.sample_ids.at(current_seg_idx);
      // For now, we do not count unphased variants as a part of this,
      // so we verify at least one of the het-pair is a "1",
      // (i.e., "-7") is treated as "0".
      // More work is needed to verify inclusion of unphased variants helps at all
      if (genotype.at(sample) != genotype.at(target_id) &&
          (genotype.at(sample) == 1 || genotype.at(target_id) == 1)) {
        path.het_sites.push_back(het_sites_processed);
      }
    }
  }
  het_sites_processed++;
}

void ThreadsLowMem::date_segments() {
  if (het_sites_processed != num_sites) {
    throw std::runtime_error(
        "Can't date segments, not all sites have been parsed for heterozygosity.");
  }
  // for (int i = 1; i < num_samples; i++) {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    if (segment_indices.at(target_id) != paths.at(target_id).size() - 1) {
      std::string prompt =
          "incomplete path at sample " + std::to_string(target_id) + ", processed ";
      prompt += std::to_string(segment_indices.at(target_id) + 1) + " segments, expected ";
      prompt += std::to_string(paths.at(target_id).size());
      throw std::runtime_error(prompt);
    }
  }

  // I get segfault on this?

  // FIXME Alex/Arni - in progress, but #def couts for now?
  // for (int i = 1; i < num_samples; i++) {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    // cout << "dating sample " << i << endl;
    // compute length in bp/cM and num_hets for each segment
    ViterbiPath& path = paths.at(target_id);
    ViterbiPath new_path(target_id);
    std::size_t n_segs = path.segment_starts.size();
    for (std::size_t k = 0; k < n_segs; k++) {
      int sample_id = path.sample_ids.at(k);
      int segment_start = path.segment_starts.at(k);
      int segment_end = k < n_segs - 1 ? path.segment_starts.at(k + 1) : num_sites - 1;
      // cout << "dating segment [" << segment_start << ", " << segment_end << ") no. " << k << "
      // out of " << n_segs - 1 << endl;
      if (segment_end == segment_start) {
        continue;
      }

      // This is inefficient but probably not that bad
      std::vector<int> segment_hets;
      for (int h : path.het_sites) {
        // FIXME Alex/Arnie review parenthesis
        if (((segment_start <= h) && (h < segment_end)) ||
            ((h == num_sites - 1) && (segment_end == num_sites - 1))) {
          segment_hets.push_back(h);
        }
        else if (h >= segment_end) {
          break;
        }
      }

      if ((target_id < n_hmm_samples) && (static_cast<int>(segment_hets.size()) > hmm_min_sites)) {
        std::vector<bool> het_hom_sites(segment_end - segment_start, false);
        for (int h : segment_hets) {
          het_hom_sites[h - segment_start] = true;
        }
        // Here we use the hmm to break the big segment up into smaller segments
        std::vector<int> breakpoints = psmc.breakpoints(het_hom_sites, segment_start);
        // cout << "found " << breakpoints.size() << " breakpoints" << endl;
        for (std::size_t j = 0; j < breakpoints.size(); j++) {
          int breakpoint_start = breakpoints[j];
          int breakpoint_end = (j == breakpoints.size() - 1) ? segment_end : breakpoints[j + 1];
          // there may be off-by-one errors here on the last segment (but who cares?)
          double bp_size =
              physical_positions.at(breakpoint_end) - physical_positions.at(breakpoint_start);
          double cm_size =
              genetic_positions.at(breakpoint_end) - genetic_positions.at(breakpoint_start);

          // Same as above
          std::vector<int> breakpoint_hets;
          for (int h : segment_hets) {
            if (((breakpoint_start <= h) && (h < breakpoint_end)) ||
                ((h == num_sites - 1) && (breakpoint_end == num_sites - 1))) {
              breakpoint_hets.push_back(h);
            }
            else if (h >= breakpoint_end) {
              break;
            }
          }
          double height = ThreadsFastLS::date_segment(static_cast<int>(breakpoint_hets.size()),
                                                      cm_size, bp_size, mutation_rate, demography);
          new_path.append(breakpoint_start, sample_id, height, breakpoint_hets);
        }
      }
      else {
        // there are off-by-one errors here on the last segment (but who cares?)
        double bp_size = physical_positions.at(segment_end) - physical_positions.at(segment_start);
        double cm_size = genetic_positions.at(segment_end) - genetic_positions.at(segment_start);
        double height = ThreadsFastLS::date_segment(
            static_cast<int>(segment_hets.size()), cm_size, bp_size, mutation_rate, demography);
        new_path.append(segment_start, sample_id, height, segment_hets);
      }
    }
    paths.at(target_id) = new_path;
  }
  return;
}

int ThreadsLowMem::count_branches() {
  int n_branches = 0;
  // for (int i = 1; i < num_samples; i++) {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    n_branches += hmms.at(target_id).count_branches();
  }
  return n_branches;
}

void ThreadsLowMem::prune() {
  // for (int i = 1; i < num_samples; i++) {
  for (int target_id : target_ids) {
    if (target_id == 0) {
      continue;
    }
    hmms.at(target_id).prune();
  }
}

// a tuple of
// segment_starts, sample_ids, heights, het_sites
std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
           std::vector<std::vector<double>>, std::vector<std::vector<int>>>
ThreadsLowMem::serialize_paths() {
  std::vector<std::vector<int>> all_starts;
  all_starts.reserve(paths.size());
  std::vector<std::vector<int>> all_ids;
  all_ids.reserve(paths.size());
  std::vector<std::vector<double>> all_heights;
  all_heights.reserve(paths.size());
  std::vector<std::vector<int>> all_hetsites;
  all_hetsites.reserve(paths.size());
  for (const int target_id : target_ids) {
    all_starts.push_back(paths.at(target_id).segment_starts);
    all_ids.push_back(paths.at(target_id).sample_ids);
    all_heights.push_back(paths.at(target_id).heights);
    all_hetsites.push_back(paths.at(target_id).het_sites);
  }
  return std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
                    std::vector<std::vector<double>>, std::vector<std::vector<int>>>(
      all_starts, all_ids, all_heights, all_hetsites);
}