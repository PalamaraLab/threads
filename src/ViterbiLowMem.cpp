#include "ViterbiLowMem.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

TracebackNode::TracebackNode(int _sample_id, int _site, TracebackNode* _previous, double _score)
    : sample_id(_sample_id), site(_site), previous(_previous), score(_score) {
}

size_t TracebackNode::key() {
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

int ViterbiPath::size() {
  return segment_starts.size();
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

void ViterbiPath::map_positions(std::vector<int>& positions) {
  bp_starts.reserve(size());
  for (auto s : segment_starts) {
    bp_starts.push_back(positions.at(s));
  }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>
ViterbiPath::dump_data_in_range(int start, int end) {
  int n_segs = size();
  if (start == -1 && end == -1 || n_segs == 0) {
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
    if (seg_start <= tmp_start && tmp_start < seg_end ||
        tmp_start <= seg_start && seg_start < tmp_end) {
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
  for (int sample_id : sample_ids) {
    size_t key = coord_id_key(0, sample_id);
    traceback_states.emplace(key, TracebackNode(sample_id, 0, nullptr, 0.));
    current_tracebacks[sample_id] = &traceback_states.at(key);
  }
  best_score = 0;
  best_match = sample_ids.at(0);
}

// Something where -9 is "missing" and -7 is "unphased"
void ViterbiState::process_site(const std::vector<int>& genotype, double rho, double rho_c,
                                double mu, double mu_c) {
  int current_site = sites_processed; // + 1;
  // be smart here! use rho_c + best_state_copy_penalty! most often this will be enough, will save a
  // lot of trouble
  double best_new_score = best_score + std::max(rho, rho_c) + std::max(mu, mu_c);
  int best_new_match = best_match;
  double new_score;
  int observed_allele = genotype.at(target_id);
  TracebackNode* prev_best = current_tracebacks.at(best_match);
  for (int sample_id : sample_ids) {
    int allele = genotype.at(sample_id);
    double copy_penalty;
    // "-7" encodes unphased heterozygotes, because why not
    if (allele == -7 || observed_allele == -7) {
      copy_penalty = (mu_c + mu) / 2.;
    }
    else {
      copy_penalty = (allele == observed_allele) ? mu_c : mu;
    }
    if (!current_tracebacks.count(sample_id)) {
      // If we've just added new sites (this will happen vary rarely),
      // recombine from previous best state
      new_score = best_score + copy_penalty + rho;
      size_t key = coord_id_key(current_site, sample_id);
      traceback_states.emplace(key, TracebackNode(sample_id, current_site, prev_best, new_score));
      current_tracebacks[sample_id] = &traceback_states.at(key);
    }
    else {
      // Otherwise, check whether we should recombine or extend
      TracebackNode* state = current_tracebacks.at(sample_id);
      if (state->score + rho_c <= best_score + rho) {
        // If extending is cheaper, simply update the score of the current traceback
        new_score = state->score + copy_penalty + rho_c;
        state->score = new_score;
      }
      else {
        // If we recombine, add a new branch
        // D-R-Y!
        new_score = best_score + copy_penalty + rho;
        size_t key = coord_id_key(current_site, sample_id);
        traceback_states.emplace(key, TracebackNode(sample_id, current_site, prev_best, new_score));
        current_tracebacks.at(sample_id) = &traceback_states.at(key);
      }
    }
    if (new_score < best_new_score) {
      best_new_score = new_score;
      best_new_match = sample_id;
    }
  }
  best_score = best_new_score;
  best_match = best_new_match;
  sites_processed++;
}

void ViterbiState::set_samples(std::unordered_set<int> new_sample_ids) {
  std::vector<int> new_samples_vec(new_sample_ids.begin(), new_sample_ids.end());
  if (!new_sample_ids.count(best_match)) {
    new_samples_vec.push_back(best_match);
  }
  for (int sample_id : sample_ids) {
    // clean up branches we definitely won't use
    if (!new_sample_ids.count(sample_id) && sample_id != best_match) {
      current_tracebacks.erase(sample_id);
    }
  }
  sample_ids = new_samples_vec;
}

void ViterbiState::prune() {
  // is there a way to do this without doing it twice?
  std::unordered_map<size_t, TracebackNode> tmp_traceback_states;

  for (int sample_id : sample_ids) {
    TracebackNode* state = current_tracebacks.at(sample_id);
    TracebackNode* new_state = recursive_insert(tmp_traceback_states, state);
    current_tracebacks[sample_id] = new_state;
  }

  traceback_states.clear();
  // cout << "real prune... ";
  for (int sample_id : sample_ids) {
    TracebackNode* state = current_tracebacks.at(sample_id);
    TracebackNode* new_state = recursive_insert(traceback_states, state);
    current_tracebacks[sample_id] = new_state;
  }
}

// add everything above and return a key to the new address
TracebackNode* ViterbiState::recursive_insert(std::unordered_map<size_t, TracebackNode>& state_map,
                                              TracebackNode* state) {
  if (state == nullptr) {
    return nullptr;
  }
  size_t key = state->key();
  if (!state_map.count(key)) {
    TracebackNode* parent = state->previous;
    TracebackNode* new_parent = recursive_insert(state_map, parent);
    state_map.emplace(key, TracebackNode(state->sample_id, state->site, new_parent, state->score));
  }
  return &state_map.at(key);
}

int ViterbiState::count_branches() {
  return traceback_states.size();
}

ViterbiPath ViterbiState::traceback() {
  ViterbiPath path(target_id);
  path.score = best_score;
  TracebackNode* state = current_tracebacks.at(best_match);
  while (state != nullptr) {
    int match_id = state->sample_id;
    int seg_start = state->site;
    path.append(seg_start, match_id);
    state = state->previous;
  }
  path.reverse();
  return path;
}
