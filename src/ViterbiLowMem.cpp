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

#include "ViterbiLowMem.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

const int ALLELE_UNPHASED_HET = -7;

inline std::size_t coord_id_key(int i, int j) {
  return (static_cast<std::size_t>(static_cast<uint32_t>(i)) << 32) |
         static_cast<std::size_t>(static_cast<uint32_t>(j));
}

} // namespace

TracebackNode::TracebackNode(int _sample_id, int _site, TracebackNode* _previous, double _score)
    : sample_id(_sample_id), site(_site), previous(_previous), score(_score) {
}

std::size_t TracebackNode::key() const {
  return coord_id_key(site, sample_id);
}

ViterbiPath::ViterbiPath(int _target_id) : target_id(_target_id) {
}

ViterbiPath::ViterbiPath(int _target_id, std::vector<int> _segment_starts,
                         std::vector<int> _sample_ids, std::vector<double> _heights,
                         std::vector<int> _het_sites)
    : target_id(_target_id), segment_starts(_segment_starts), sample_ids(_sample_ids),
      heights(_heights), het_sites(_het_sites) {
}

void ViterbiPath::reverse() {
  std::reverse(segment_starts.begin(), segment_starts.end());
  std::reverse(sample_ids.begin(), sample_ids.end());
}

int ViterbiPath::size() const {
  return static_cast<int>(segment_starts.size());
}

void ViterbiPath::append(int segment_start, int sample_id) {
  segment_starts.push_back(segment_start);
  sample_ids.push_back(sample_id);
}

void ViterbiPath::append(int segment_start, int sample_id, double height,
                         std::vector<int>& new_het_sites) {
  if (segment_starts.size() > 0 && segment_start <= segment_starts.back()) {
    throw std::runtime_error("Illegal path append of segment start " +
                             std::to_string(segment_start) + " after segment " +
                             std::to_string(segment_starts.back()));
  }
  if (het_sites.size() > 0 && segment_start <= het_sites.back()) {
    throw std::runtime_error("Illegal path append of segment start " +
                             std::to_string(segment_start) + " after het-site " +
                             std::to_string(het_sites.back()));
  }
  segment_starts.push_back(segment_start);
  sample_ids.push_back(sample_id);
  heights.push_back(height);
  het_sites.insert(het_sites.end(), new_het_sites.begin(), new_het_sites.end());
}

void ViterbiPath::map_positions(const std::vector<int>& positions) {
  bp_starts.reserve(size());
  for (auto s : segment_starts) {
    bp_starts.push_back(positions.at(s));
  }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>
ViterbiPath::dump_data_in_range(int start, int end) {
  int n_segs = size();
  if (((start == -1) && (end == -1)) || (n_segs == 0)) {
    return std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>(
        bp_starts, sample_ids, heights);
  }
  // these are in base-pair units
  int tmp_end = end == -1 ? std::numeric_limits<int>::max() : end;
  int tmp_start = std::max(start, bp_starts.at(0));

  std::vector<int> out_starts;
  std::vector<int> out_ids;
  std::vector<double> out_heights;
  for (int i = 0; i < n_segs; i++) {
    int seg_start = bp_starts.at(i);
    int seg_end = i < n_segs - 1 ? bp_starts.at(i + 1) : tmp_end;
    if (((seg_start <= tmp_start) && (tmp_start < seg_end)) ||
        ((tmp_start <= seg_start) && (seg_start < tmp_end))) {
      out_starts.push_back(std::max(tmp_start, seg_start));
      out_ids.push_back(sample_ids.at(i));
      out_heights.push_back(heights.at(i));
    }
    else if (tmp_end <= seg_start) {
      break;
    }
  }
  return std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>(
      out_starts, out_ids, out_heights);
}

ViterbiState::ViterbiState(int _target_id, std::vector<int> _sample_ids)
    : target_id(_target_id), sample_ids(_sample_ids) {
  // init current_tracebacks
  if (sample_ids.size() == 0) {
    throw std::runtime_error("found no samples for ViterbiState object for sample " +
                             std::to_string(target_id));
  }
  current_traceback_ptrs.reserve(sample_ids.size());
  for (int sample_id : sample_ids) {
    current_traceback_ptrs.push_back(alloc_node(sample_id, 0, nullptr, 0.));
  }
  best_score = 0;
  best_match = sample_ids.at(0);
  best_match_idx = 0;
}

void ViterbiState::process_site(const int* genotype, double rho, double rho_c,
                                double mu, double mu_c) {
  const int current_site = sites_processed;
  double best_new_score = best_score + std::max(rho, rho_c) + std::max(mu, mu_c);
  int best_new_match = best_match;
  int best_new_match_idx = best_match_idx;
  const int observed_allele = genotype[target_id];
  const double recomb_threshold = best_score + rho;
  const double unphased_penalty = (mu_c + mu) * 0.5;
  const bool observed_is_unphased = (observed_allele == ALLELE_UNPHASED_HET);

  TracebackNode* prev_best = current_traceback_ptrs[best_match_idx];

  const int n_samples = static_cast<int>(sample_ids.size());
  for (int idx = 0; idx < n_samples; ++idx) {
    const int sample_id = sample_ids[idx];
    const int allele = genotype[sample_id];
    double copy_penalty;
    if (observed_is_unphased || (allele == ALLELE_UNPHASED_HET)) {
      copy_penalty = unphased_penalty;
    }
    else {
      copy_penalty = (allele == observed_allele) ? mu_c : mu;
    }

    double new_score;
    TracebackNode* state = current_traceback_ptrs[idx];
    if (state == nullptr) {
      // Newly added sample (happens rarely, after set_samples)
      new_score = recomb_threshold + copy_penalty;
      current_traceback_ptrs[idx] = alloc_node(sample_id, current_site, prev_best, new_score);
    }
    else {
      if (state->score + rho_c <= recomb_threshold) {
        // Extend: cheaper than recombining
        new_score = state->score + copy_penalty + rho_c;
        state->score = new_score;
      }
      else {
        // Recombine: add a new branch
        new_score = recomb_threshold + copy_penalty;
        current_traceback_ptrs[idx] = alloc_node(sample_id, current_site, prev_best, new_score);
      }
    }
    if (new_score < best_new_score) {
      best_new_score = new_score;
      best_new_match = sample_id;
      best_new_match_idx = idx;
    }
  }
  best_score = best_new_score;
  best_match = best_new_match;
  best_match_idx = best_new_match_idx;
  sites_processed++;
}

void ViterbiState::set_samples(std::unordered_set<int> new_sample_ids) {
  // Build old sample_id → ptr map from current parallel vectors
  std::unordered_map<int, TracebackNode*> old_ptrs;
  old_ptrs.reserve(sample_ids.size());
  for (std::size_t i = 0; i < sample_ids.size(); ++i) {
    old_ptrs[sample_ids[i]] = current_traceback_ptrs[i];
  }

  std::vector<int> new_samples_vec(new_sample_ids.begin(), new_sample_ids.end());
  if (!new_sample_ids.count(best_match)) {
    new_samples_vec.push_back(best_match);
  }
  sample_ids = new_samples_vec;

  // Rebuild parallel pointer vector to match new sample_ids ordering
  current_traceback_ptrs.clear();
  current_traceback_ptrs.reserve(sample_ids.size());
  for (std::size_t i = 0; i < sample_ids.size(); ++i) {
    int sample_id = sample_ids[i];
    auto it = old_ptrs.find(sample_id);
    current_traceback_ptrs.push_back(it != old_ptrs.end() ? it->second : nullptr);
    if (sample_id == best_match) {
      best_match_idx = static_cast<int>(i);
    }
  }
}

void ViterbiState::prune() {
  std::deque<TracebackNode> new_nodes;
  std::unordered_map<std::size_t, TracebackNode*> key_to_ptr;

  // Recursively copy only reachable nodes into the new deque
  auto copy_node = [&](auto& self, TracebackNode* state) -> TracebackNode* {
    if (state == nullptr) return nullptr;
    std::size_t key = state->key();
    auto it = key_to_ptr.find(key);
    if (it != key_to_ptr.end()) return it->second;
    TracebackNode* new_parent = self(self, state->previous);
    new_nodes.emplace_back(state->sample_id, state->site, new_parent, state->score);
    TracebackNode* ptr = &new_nodes.back();
    key_to_ptr[key] = ptr;
    return ptr;
  };

  for (std::size_t idx = 0; idx < sample_ids.size(); ++idx) {
    current_traceback_ptrs[idx] = copy_node(copy_node, current_traceback_ptrs[idx]);
  }

  traceback_nodes = std::move(new_nodes);
}

// add everything above and return a key to the new address
TracebackNode*
ViterbiState::recursive_insert(std::unordered_map<std::size_t, TracebackNode>& state_map,
                               TracebackNode* state) {
  if (state == nullptr) {
    return nullptr;
  }
  std::size_t key = state->key();
  if (!state_map.count(key)) {
    TracebackNode* parent = state->previous;
    TracebackNode* new_parent = recursive_insert(state_map, parent);
    state_map.emplace(key, TracebackNode(state->sample_id, state->site, new_parent, state->score));
  }
  return &state_map.at(key);
}

TracebackNode* ViterbiState::alloc_node(int sample_id, int site, TracebackNode* previous, double score) {
  traceback_nodes.emplace_back(sample_id, site, previous, score);
  return &traceback_nodes.back();
}

int ViterbiState::count_branches() const {
  return static_cast<int>(traceback_nodes.size());
}

ViterbiPath ViterbiState::traceback() {
  ViterbiPath path(target_id);
  path.score = best_score;
  TracebackNode* state = current_traceback_ptrs[best_match_idx];
  while (state != nullptr) {
    int match_id = state->sample_id;
    int seg_start = state->site;
    path.append(seg_start, match_id);
    state = state->previous;
  }
  path.reverse();
  return path;
}
