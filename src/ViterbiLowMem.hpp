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

#ifndef THREADS_ARG_VITERBI_LOW_MEM_HPP
#define THREADS_ARG_VITERBI_LOW_MEM_HPP

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class TracebackNode {
public:
  TracebackNode(int _sample_id, int _site, TracebackNode* _previous, double _score);
  size_t key();

public:
  int sample_id = 0;
  int site = 0;
  TracebackNode* previous = nullptr;
  double score = 0.0;
};

class ViterbiPath {
public:
  ViterbiPath(int _target_id);
  ViterbiPath(int _target_id, std::vector<int> _segment_starts, std::vector<int> _sample_ids,
              std::vector<double> _heights, std::vector<int> _het_sites);
  void append(int segment_start, int sample_id);
  void append(int segment_start, int sample_id, double height, std::vector<int>& new_het_sites);
  std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>
  dump_data_in_range(int start_idx, int end_idx);
  void reverse();
  int size() const;
  void map_positions(std::vector<int>& positions);

public:
  double score = 0.0;
  int target_id = 0;
  std::vector<int> bp_starts;
  std::vector<int> segment_starts;
  std::vector<int> sample_ids;
  std::vector<double> heights;
  std::vector<int> het_sites;
};

class ViterbiState {
public:
  ViterbiState(int _target_id, std::vector<int> _sample_ids);

  void process_site(const std::vector<int>& genotype, double rho, double rho_c, double _mu,
                    double _mu_c);
  void set_samples(std::unordered_set<int> new_sample_ids);
  int count_branches() const;
  void prune();
  ViterbiPath traceback();

private:
  std::unordered_map<size_t, TracebackNode> traceback_states;
  TracebackNode* recursive_insert(std::unordered_map<size_t, TracebackNode>& state_map,
                                  TracebackNode* state);

public:
  int target_id = 0;
  int best_match = -1;
  double best_score = 0.0;
  int sites_processed = 0;
  double mutation_penalty = 0.0;
  std::vector<int> sample_ids;
  std::vector<double> sample_scores;
  std::unordered_map<int, TracebackNode*> current_tracebacks;
};

#endif // THREADS_ARG_VITERBI_LOW_MEM_HPP
