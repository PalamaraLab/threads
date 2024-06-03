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

#ifndef THREADS_ARG_THREADS_LOW_MEM_HPP
#define THREADS_ARG_THREADS_LOW_MEM_HPP

#include "Demography.hpp"
#include "Matcher.hpp"
#include "ThreadsFastLS.hpp"
#include "ViterbiLowMem.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ThreadsLowMem {
public:
  ThreadsLowMem(const std::vector<int> _target_ids, const std::vector<double>& _physical_positions,
                const std::vector<double>& _genetic_positions, std::vector<double> ne,
                std::vector<double> ne_times, double _mutation_rate, bool _sparse);

  // Algorithm outline

  // 1. process all sites for the PBWT (done by the Matcher)

  // 2a. initialize hmms
  void initialize_viterbi(std::vector<std::vector<std::unordered_set<int>>>& match_ids,
                          const std::vector<double>& cm_positions);
  // 2b. process all sites for the hmms
  void process_site_viterbi(const std::vector<int>& genotype);
  // 2c. prune branches at regular intervals (i.e. when there's a lot of them, figure this out soon)
  void prune();
  // 2d. traceback all the hmms to get viterbi paths
  void traceback();

  // 3a. add het sites
  void process_site_hets(const std::vector<int>& genotype);
  // 3b. date all segments
  void date_segments();

  // 4. save output (done on Python side)

  int count_branches() const;

  // Make pickle-able path output
  std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
             std::vector<std::vector<double>>, std::vector<std::vector<int>>>
  serialize_paths();

public:
  // This object will only run the HMM for these ids
  std::vector<int> target_ids;
  std::unordered_map<int, double> expected_branch_lengths;
  double mean_bp_size = 0.0;
  std::unordered_map<int, int> segment_indices;
  std::unordered_map<int, ViterbiPath> paths;
  int num_samples = 0;
  int num_sites = 0;
  double mutation_rate = 0.0;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  std::vector<double> bp_boundaries;
  std::vector<double> cm_boundaries;
  bool sparse = false;

private:
  Demography demography;

  // 2. HMM quantites
  int hmm_sites_processed = 0;
  std::unordered_map<int, ViterbiState> hmms;
  int match_group_idx = 0;
  std::vector<MatchGroup> match_groups;

  // 3. Path segment and path dating quantites
  HMM psmc;
  int het_sites_processed = 0;
  int n_hmm_samples = 100;
  int hmm_min_sites = 10;
};

#endif // THREADS_ARG_THREADS_LOW_MEM_HPP
