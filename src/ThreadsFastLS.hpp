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

#ifndef THREADS_ARG_THREADS_HPP
#define THREADS_ARG_THREADS_HPP

#include "HMM.hpp"
#include "State.hpp"

#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

struct ImputationSegment {
  int seg_start = 0;
  std::vector<int> ids;

  // these are at the *start* of the segment
  std::vector<double> weights;
  std::vector<double> ages;
};

/// This class contains everything for Threads-fastLS
class ThreadsFastLS {
public:
  // Constructors and utils
  ThreadsFastLS(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
                double _mutation_rate, double ne, bool _sparse_sites, int _n_prune, bool _use_hmm,
                int _burn_in_left, int _burn_in_right)
      : ThreadsFastLS(_physical_positions, _genetic_positions, _mutation_rate,
                      std::vector<double>{ne}, std::vector<double>{0.0}, _sparse_sites, _n_prune,
                      _use_hmm, _burn_in_left, _burn_in_right){};
  ThreadsFastLS(std::vector<double> _physical_positions, std::vector<double> _genetic_positions,
                double _mutation_rate, std::vector<double> ne, std::vector<double> ne_times,
                bool _sparse_sites, int _n_prune, bool _use_hmm, int _burn_in_left,
                int _burn_in_right);

  /// Find the next insert position
  /// @param t Pointer to node below the sequence being inserted at site i
  /// @param g Allele at site i+1
  /// @param i The site
  /// @return Node* Pointer to node below the sequence being inserted at site i+1
  Node* extend_node(Node* node, bool genotype, const int i);

  /// Determine whether state can be extended through panel by appending genotype 'g'
  /// @param s State at site i
  /// @param t_next Node at site i+1
  /// @param g Candidate genotype for s at i+1
  /// @param i The site index
  /// @return bool true if extensible
  bool extensible_by(State& s, const Node* t_next, const bool g, const int i);

  std::pair<bool, bool> pair_extensible_by(StatePair& p, const bool g, const int i);

  /// @param id1
  /// @param id2
  /// @param start inclusive!
  /// @param end exclusive!
  /// @return bool whether sequences match on the interval
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);

  /// Assuming input genotypes match sample_id on [segment_start, segment_end), how much can we
  /// extend the region in either direction with out hitting a mismatch
  std::pair<int, int> overflow_region(const std::vector<bool>& genotypes, const int sample_id,
                                      const int segment_start, const int segment_end);

  /// Fetch het-hom status for id1 and id2 in the region specified by site indices
  std::vector<bool> fetch_het_hom_sites(const int id1, const int id2, const int start,
                                        const int end);

  std::vector<int> het_sites_from_thread(const int focal_ID, std::vector<int> bp_starts,
                                         std::vector<std::vector<int>> target_IDs);

  static std::tuple<std::vector<double>, std::vector<double>>
  site_sizes(std::vector<double> positions);

  // More attributes
  std::vector<double> trimmed_positions() const;

  /// Insert and assign generic ID
  void insert(const std::vector<bool>& genotype);

  /// Insert a new sequence into the dynamic panel and assign specific ID
  void insert(const int ID, const std::vector<bool>& genotype);

  /// Deletes sequence ID from the dynamic panel. This moves the last sequence in the panel
  /// to the position ID held. See alg 5 from d-PBWT paper.
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

  std::vector<ImputationSegment> impute(std::vector<bool>& genotype, int neighborhood_size);

  /// Run Li-Stephens on input haplotype *without* inserting into the dynamic panel.
  /// See also Algorithm 4 of Lunter (2018), Bioinformatics.
  /// For imputation we use the IMPUTE/Beagle mutation penalties
  std::pair<TracebackState*, Node*> fastLS(const std::vector<bool>& genotype,
                                           bool imputation = false);

  /// This is a basic traceback that samples a random haplotype per from all matches per segment and
  /// also stores the lowest-numbered match (for compression)
  std::vector<std::tuple<int, std::vector<int>>> traceback(TracebackState* tb, Node* match,
                                                           bool return_all = false);

  /// Similar to normal traceback, but picks the up-to neighborhood_size best matches and stores
  /// their overlap with the input sequence returns a list of a tuple of lists :-P I.e., a list of
  /// segments, and each segment is a tuple (sample_IDs, starts, ends) all of equal length <=
  /// neighborhood_size
  std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>>
  traceback_impute(std::vector<bool>& genotypes, TracebackState* tb, Node* match,
                   int neighborhood_size);

  /// Run Li-Stephens on input diplotype *without* inserting into the dynamic panel.
  /// See also Algorithm 4 of Lunter (2018), Bioinformatics.
  /// Warning: The code here is more verbose than it has to be
  std::array<std::pair<TracebackState*, Node*>, 2> fastLS_diploid(const std::vector<int>& genotype);

  std::array<std::vector<std::tuple<int, std::vector<int>>>, 2>
  diploid_ls(std::vector<int> unphased_genotypes);

  /// This gives the IMPUTE5 (and Beagle) recombination penalties
  /// @return tuple of penalty vectors
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties_impute5();

  /// This gives the *sparse* recombination penalties
  /// @return tuple of penalty vectors
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties();
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties();

  /// This gives the correct, dense, recombination penalties
  /// @return tuple of penalty vectors
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties_correct();

  // TMRCA estimation
  /// Date the segment based on length and n_mismatches using maximum likelihood. (No
  /// demography)
  /// @param id1
  /// @param id2
  /// @param start inclusive
  /// @param end exclusive
  /// @return double
  double date_segment(const int num_het_sites, const int start, const int end);

  static double date_segment(int num_het_sites, double cm_length, double bp_length,
                             double mutation_rate, Demography& demography);
  static double date_segment_sparse(double cm_length, Demography& demography);

  /// For debugging: print the sample-IDs of the arrayified panel
  void print_sorting() const;

public:
  int n_prune = 0;
  int num_sites = 0;
  int num_samples = 0;
  double mutation_rate = 0.0;
  int burn_in_left = 0;
  int burn_in_right = 0;
  double threading_start = 0.0;
  double threading_end = 0.0;

  bool sparse_sites = false; ///< This determines the segment dating formula
  bool use_hmm = false;      ///< Whether to internally break up big segments

  std::unordered_map<int, int> ID_map;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  std::vector<double> bp_boundaries;
  std::vector<double> cm_boundaries;
  Demography demography;
  HMM* hmm = nullptr;

  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  std::vector<std::unique_ptr<TracebackState>> traceback_states;

  // Burn-in quantities
  int trim_pos_start_idx = 0;
  int trim_pos_end_idx = 0;

  // TODO set random seed at runtime (ticket #20)
  std::mt19937 rng{1234};
};

#endif // THREADS_ARG_THREADS_HPP
