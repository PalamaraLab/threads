#include "TGEN.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;

TGEN::TGEN(std::vector<int> _positions, std::vector<std::vector<int>> _bp_starts,
           std::vector<std::vector<int>> _target_IDs, std::vector<std::vector<int>> _het_sites)
    : positions(_positions), bp_starts(_bp_starts), target_IDs(_target_IDs), het_sites(_het_sites) {
  positions.push_back(std::numeric_limits<int>::max());
  // position-to-site map
  for (int i = 0; i < positions.size(); i++) {
    pos_idx_map[positions[i]] = i;
  }

  // set reference genome here
  reference_genome = Eigen::VectorXi::Zero(positions.size());
  reference_genome_vec = std::vector<bool>(positions.size(), false);
  for (int h : het_sites[0]) {
    reference_genome[pos_idx_map.at(h)] = 1;
    reference_genome_vec[pos_idx_map.at(h)] = true;
  }

  interval_sets.reserve(bp_starts.size());
  // Initialize interval maps for each sample. This can get too slow
  for (int i = 0; i < bp_starts.size(); i++) {
    SegmentSet iset;
    int n_segs = bp_starts[i].size();
    std::vector<int>& sample_hets = het_sites[i];
    int n_hets = sample_hets.size();
    int het_site_idx = 0;
    int pos_idx = 0;
    std::vector<int> seg_hets;
    for (int j = 0; j < n_segs; j++) {
      int seg_start = bp_starts[i][j];
      int seg_end = j == n_segs - 1 ? positions.back() : bp_starts[i][j + 1];

      // Binary search through positions to get segment boundaries
      int seg_start_pos = *std::lower_bound(positions.begin(), positions.end(), seg_start);
      int seg_end_pos = j == n_segs - 1 ? positions.back()
                                        : *std::lower_bound(positions.begin(), positions.end(),
                                                            bp_starts[i][j + 1]);

      while (het_site_idx < n_hets && seg_start <= sample_hets[het_site_idx] &&
             sample_hets[het_site_idx] < seg_end) {
        seg_hets.push_back(sample_hets[het_site_idx]);
        het_site_idx++;
      }
      iset += TgenSegment(seg_start_pos, seg_end_pos, seg_hets, target_IDs[i][j]);
    }
    interval_sets.push_back(iset);
  }

  // initialize cached genotypes
  cached_genotypes_map[0] = -1;
}

// Eigen-based query
Eigen::MatrixXi& TGEN::query(const int bp_from, const int bp_to, const std::vector<int>& samples) {
  clear_cache();

  // Find number of expected sites
  int start_pos = *std::lower_bound(positions.begin(), positions.end(), bp_from);
  int end_pos = *std::upper_bound(positions.begin(), positions.end(), bp_to);
  int idx_offset = pos_idx_map[start_pos];
  genotype_cache.resize(samples.size(), pos_idx_map[end_pos] - idx_offset);

  TgenSegment range(start_pos, end_pos);

  for (int i = 0; i < samples.size(); i++) {
    cached_genotypes_map[samples[i]] = i;
    if (samples[i] == 0) {
      auto insert_range =
          Eigen::seq(0, pos_idx_map[end_pos] - idx_offset - 1); // eigen seq is inclusive
      auto copy_range = Eigen::seq(idx_offset, pos_idx_map[end_pos] - 1);
      genotype_cache(i, insert_range) = reference_genome(
          copy_range); //.eval(); //WARNING need .eval() here (or do we? I don't think we do)
    }
    else {
      SegmentSet& segments(interval_sets[samples[i]]);

      // Initialize the queue
      std::queue<TgenSegment> seg_queue;
      auto eqr = segments.equal_range(range);
      for (SegmentSet::const_iterator iter = eqr.first; iter != eqr.second; iter++) {
        seg_queue.push(*iter & range);
      }

      // Process everything in the queue
      while (!seg_queue.empty()) {
        TgenSegment& segment = seg_queue.front();
        if (cached_genotypes_map.find(segment.target) != cached_genotypes_map.end()) {
          // We've reached somewhere along the tree where we can copy from
          int seg_start_idx = pos_idx_map[segment.lower()];
          int seg_end_idx = pos_idx_map[segment.upper()];
          auto insert_range = Eigen::seq(
              seg_start_idx - idx_offset, seg_end_idx - idx_offset - 1); // eigen seq is inclusive
          if (segment.target == 0) {
            auto copy_range = Eigen::seq(seg_start_idx, seg_end_idx - 1);
            // We've reached the root of the tree and copy from the "reference" genome
            genotype_cache(i, insert_range) = reference_genome(copy_range);
          }
          else {
            // We've found a cached genotype to copy from

            genotype_cache(i, insert_range) =
                genotype_cache(cached_genotypes_map[segment.target], insert_range);
          }

          // We then flip all the het sites
          for (int h : segment.het_sites) {
            genotype_cache(i, pos_idx_map[h] - idx_offset) =
                1 - genotype_cache(i, pos_idx_map[h] - idx_offset);
          }
        }
        else {
          // We've not yet reached somewhere to copy from, so we keep traversing
          auto new_eqr = interval_sets[segment.target].equal_range(segment);
          for (SegmentSet::const_iterator iter = new_eqr.first; iter != new_eqr.second; iter++) {
            seg_queue.push(*iter & segment);
          }
        }
        seg_queue.pop();
      }
    }
  }
  return genotype_cache;
}

// std::vector-based query
// Warning: This makes a copy when returned through the python interface
std::vector<std::vector<bool>>& TGEN::query2(const int bp_from, const int bp_to,
                                             const std::vector<int>& samples) {
  genotypes.clear();

  // Find number of expected sites
  int start_pos = *std::lower_bound(positions.begin(), positions.end(), bp_from);
  int end_pos = *std::upper_bound(positions.begin(), positions.end(), bp_to);
  int idx_offset = pos_idx_map[start_pos];

  int n_samples = samples.size();
  int n_sites = pos_idx_map[end_pos] - idx_offset;
  genotypes.reserve(samples.size());
  for (int i = 0; i < samples.size(); i++) {
    genotypes.push_back(std::vector<bool>(n_sites));
  }
  TgenSegment range(start_pos, end_pos);

  for (int i = 0; i < samples.size(); i++) {
    std::vector<bool>& current_gt = genotypes.at(i);
    cached_genotypes_map[samples[i]] = i;
    if (samples[i] == 0) {
      std::copy(reference_genome_vec.begin() + idx_offset,
                reference_genome_vec.begin() + pos_idx_map[end_pos], current_gt.begin());
    }
    else {
      SegmentSet& segments(interval_sets[samples[i]]);

      // Initialize the queue
      std::queue<TgenSegment> seg_queue;
      auto eqr = segments.equal_range(range);
      for (SegmentSet::const_iterator iter = eqr.first; iter != eqr.second; iter++) {
        seg_queue.push(*iter & range);
      }

      // Process everything in the queue
      while (!seg_queue.empty()) {
        TgenSegment& segment = seg_queue.front();
        if (cached_genotypes_map.find(segment.target) != cached_genotypes_map.end()) {
          // We've reached somewhere along the tree where we can copy from
          int seg_start_idx = pos_idx_map[segment.lower()];
          int seg_end_idx = pos_idx_map[segment.upper()];
          if (segment.target == 0) {
            // We've reached the root of the tree and copy from the "reference" genome
            std::copy(reference_genome_vec.begin() + seg_start_idx,
                      reference_genome_vec.begin() + seg_end_idx,
                      current_gt.begin() + seg_start_idx - idx_offset);
          }
          else {
            std::vector<bool>& target_gt = genotypes.at(cached_genotypes_map[segment.target]);
            std::copy(target_gt.begin() + seg_start_idx - idx_offset,
                      target_gt.begin() + seg_end_idx - idx_offset,
                      current_gt.begin() + seg_start_idx - idx_offset);
            // We've found a cached genotype to copy from
          }
          // We then flip all the het sites
          // then we can also delay the .eval() until the end, right?
          for (int h : segment.het_sites) {
            current_gt.at(pos_idx_map[h] - idx_offset) =
                !current_gt.at(pos_idx_map[h] - idx_offset);
          }
        }
        else {
          // We've not yet reached somewhere to copy from, so we keep traversing
          auto new_eqr = interval_sets[segment.target].equal_range(segment);
          for (SegmentSet::const_iterator iter = new_eqr.first; iter != new_eqr.second; iter++) {
            seg_queue.push(*iter & segment);
          }
        }
        seg_queue.pop();
      }
    }
  }
  return genotypes;
}

void TGEN::clear_cache() {
  cached_genotypes_map.clear();
  genotype_cache.resize(0, 0);
  cached_genotypes_map[0] = -1;
}
