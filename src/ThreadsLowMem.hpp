#ifndef THREADS_INFER_THREADS_LOW_MEM_HPP
#define THREADS_INFER_THREADS_LOW_MEM_HPP

#include "Matcher.hpp"
#include "Threads.hpp"
#include "ViterbiLowMem.hpp"
#include "Demography.hpp"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ThreadsLowMem {

private:
  Demography demography;

  // 2. HMM quantites
  int hmm_sites_processed;
  std::unordered_map<int, ViterbiState> hmms;
  int match_group_idx;
  std::vector<MatchGroup> match_groups;

  // 3. Path segment and path dating quantites
  HMM psmc;
  // std::vector<double> cm_sizes;
  // std::vector<double> bp_sizes;
  int het_sites_processed;
  int n_hmm_samples = 100;
  int hmm_min_sites = 10;

public:
  // This object will only run the HMM for these ids
  std::vector<int> target_ids;
  std::unordered_map<int, double> expected_branch_lengths;
  double mean_bp_size;
  std::unordered_map<int, int> segment_indices;
  std::unordered_map<int, ViterbiPath> paths;
  int num_samples;
  int num_sites;
  double mutation_rate;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  std::vector<double> bp_boundaries;
  std::vector<double> cm_boundaries;
  bool sparse;

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

  int count_branches();

  // Make pickle-able path output
  std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>,
             std::vector<std::vector<double>>, std::vector<std::vector<int>>>
  serialize_paths();
};

#endif // THREADS_INFER_THREADS_LOW_MEM_HPP
