#ifndef THREADS_INFER_THREADS_HPP
#define THREADS_INFER_THREADS_HPP

#include "State.hpp"
#include "HMM.hpp"

#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

struct ImputationSegment {
  int seg_start;
  std::vector<int> ids;
  // these are at the *start* of the segment
  std::vector<double> weights;
  std::vector<double> ages;
};

// This class contains everything for Threads-fastLS
class Threads {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;  // Ha, ha

  std::vector<std::unique_ptr<TracebackState>> traceback_states;

  // Burn-in quantities
  int trim_pos_start_idx;
  int trim_pos_end_idx;

  std::mt19937 rng;

  inline size_t pair_key(int i, int j) {
    return (size_t) i << 32 | (unsigned int) j;
  }
  Node* extend_node(Node* node, bool genotype, const int i);
  bool extensible_by(State& s, const Node* t_next, const bool g, const int i);
  std::pair<bool, bool> pair_extensible_by(StatePair& p, const bool g, const int i);
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);
  std::pair<int, int> overflow_region(const std::vector<bool>& genotypes, const int sample_id,
                                      const int segment_start, const int segment_end);
  std::vector<bool> fetch_het_hom_sites(const int id1, const int id2, const int start, const int end);
  std::vector<int> het_sites_from_thread(const int focal_ID, std::vector<int> bp_starts,
                                         std::vector<std::vector<int>> target_IDs);

public:
  int n_prune;
  int num_sites;
  int num_samples;
  double mutation_rate;
  int burn_in_left;
  int burn_in_right;
  double threading_start;
  double threading_end;
  // This determines the segment dating formula
  bool sparse_sites;
  // Whether to internally break up big segments
  bool use_hmm;
  std::unordered_map<int, int> ID_map;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  std::vector<double> bp_boundaries;
  std::vector<double> cm_boundaries;
  Demography demography;
  HMM* hmm;

  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;

  // Constructors and utils
  Threads(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
          double _mutation_rate, double ne, bool _sparse_sites, int _n_prune, bool _use_hmm,
          int _burn_in_left, int _burn_in_right)
      : Threads(_physical_positions, _genetic_positions, _mutation_rate, std::vector<double>{ne},
                std::vector<double>{0.0}, _sparse_sites, _n_prune, _use_hmm, _burn_in_left,
                _burn_in_right){};
  Threads(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
          double _mutation_rate, std::vector<double> ne, std::vector<double> ne_times,
          bool _sparse_sites, int _n_prune, bool _use_hmm, int _burn_in_left, int _burn_in_right);

  static std::tuple<std::vector<double>, std::vector<double>>
  site_sizes(std::vector<double> positions);

  // More attributes
  std::vector<double> trimmed_positions();

  // Insertion/deletion
  // Insert and assign generic ID
  void insert(const std::vector<bool>& genotype);
  // Insert and assign specific ID
  void insert(const int ID, const std::vector<bool>& genotype);
  // Remove sample from panel
  void remove(int ID);

  // HMM
  void delete_hmm();

  // Algorithms
  std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
  thread(const std::vector<bool>& genotype);
  std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
  thread(const int new_sample_ID, const std::vector<bool>& genotype);
  std::tuple<std::vector<int>, std::vector<std::vector<int>>, std::vector<double>, std::vector<int>>
  remove_burn_in(std::vector<int>& bp_starts, std::vector<std::vector<int>>& target_IDs,
                 std::vector<double>& segment_ages, std::vector<int>& het_sites);
  std::vector<std::tuple<int, std::vector<int>>> threads_ls(const std::vector<bool>& genotype);

  std::vector<ImputationSegment> impute(std::vector<bool>& genotype, int L);
  std::pair<TracebackState*, Node*> fastLS(const std::vector<bool>& genotype,
                                           bool imputation = false);
  std::vector<std::tuple<int, std::vector<int>>> traceback(TracebackState* tb, Node* match,
                                                           bool return_all = false);
  std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>>
  traceback_impute(std::vector<bool>& genotypes, TracebackState* tb, Node* match, int L);
  std::array<std::pair<TracebackState*, Node*>, 2> fastLS_diploid(const std::vector<int>& genotype);
  std::array<std::vector<std::tuple<int, std::vector<int>>>, 2>
  diploid_ls(std::vector<int> unphased_genotypes);
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties_impute5();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties();
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties_correct();

  // TMRCA estimation
  double date_segment(const int num_het_sites, const int start, const int end);
  static double date_segment(int num_het_sites, double cm_length, double bp_length,
                             double mutation_rate, Demography& demography);
  static double date_segment_sparse(int num_het_sites, double cm_length, Demography& demography);
  // Debugging
  void print_sorting();
};

#endif // THREADS_INFER_THREADS_HPP
