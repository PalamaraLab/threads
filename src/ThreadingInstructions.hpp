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

#ifndef THREADS_ARG_THREADING_INSTRUCTIONS_HPP
#define THREADS_ARG_THREADING_INSTRUCTIONS_HPP

#include "ViterbiLowMem.hpp"

#include <limits>
#include <sstream>
#include <vector>


class ThreadingInstruction {
// This is a container class for the instructions inferred for a single sequence
public:
    // Segment starts
    std::vector<int> starts;
    // Coalescence times (times to the most recent common ancestor)
    std::vector<double> tmrcas;
    // Closest cousins
    std::vector<int> targets;
    // Mismatches with closest cousin
    std::vector<int> mismatches;
    std::size_t num_segments = 0;
    std::size_t num_mismatches = 0;
    ThreadingInstruction(
        std::vector<int> _starts,
        std::vector<double> _tmrcas,
        std::vector<int> _targets,
        std::vector<int> _mismatches
    );
};

class ThreadingInstructionIterator {
// This is a container for iterating through a ThreadingInstruction object site-by-site
private:
    const ThreadingInstruction& instruction;
    const std::vector<int>& positions;
    int sites_processed = 0;
    int next_segment_start = 0;
    int next_mismatch = 0;
public:
    int num_sites = 0;
    int current_mismatch = 0;
    int current_segment = 0;
    int current_target = 0;
    double current_tmrca = 0;
    bool is_mismatch = false;
public:
    void increment_site(const int pos);
    ThreadingInstructionIterator(const ThreadingInstruction& _instruction, const std::vector<int>& _positions);
};

class ThreadingInstructions {
// This is a container class for threading instructions inferred by Threads
public:
    // Getters
    std::vector<std::vector<int>> all_starts();
    std::vector<std::vector<double>> all_tmrcas();
    std::vector<std::vector<int>> all_targets();
    std::vector<std::vector<int>> all_mismatches();

    // Constructors
    ThreadingInstructions(std::vector<ThreadingInstruction> _instructions, std::vector<int> _positions);
    ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int start, const int end, const std::vector<int>& all_positions);
    ThreadingInstructions(const std::vector<std::vector<int>>& starts,
                          const std::vector<std::vector<double>>& tmrcas,
                          const std::vector<std::vector<int>>& targets,
                          const std::vector<std::vector<int>>& mismatches,
                          const std::vector<int>& _positions, int _start, int _end);

    ThreadingInstructions sub_range(const int range_start, const int range_end) const;

    // Common operations
    std::vector<double> left_multiply(const std::vector<double>& x, bool diploid=false, bool normalize=false);
    std::vector<double> right_multiply(const std::vector<double>& x, bool diploid=false, bool normalize=false);

    // Tree-propagation multiply: O(n * n_segments + total_mismatches) per call.
    // One-time O(n * m) prepare step precomputes mismatch corrections using
    // reference-counted genotype cache (peak memory ~O(m * tree_depth)).
    // No genotype matrix retained after prepare.
    void prepare_tree_multiply();
    std::vector<double> right_multiply_tree(const std::vector<double>& x);
    std::vector<double> left_multiply_tree(const std::vector<double>& x);

    // Batch tree multiply: process k vectors in a single tree traversal.
    // Input/output are row-major flat arrays: X[row * k + col].
    // right_multiply_tree_batch: X is (num_sites, k), returns (num_samples, k)
    // left_multiply_tree_batch:  X is (num_samples, k), returns (num_sites, k)
    std::vector<double> right_multiply_tree_batch(const std::vector<double>& x_flat, int k);
    std::vector<double> left_multiply_tree_batch(const std::vector<double>& x_flat, int k);

    // Precompute and cache the dense genotype matrix (num_sites × num_samples, row-major).
    // Subsequent left_multiply/right_multiply calls use the cached matrix.
    void materialize_genotypes();

    // Precompute and cache standardized matrices for normalized multiply.
    // These store (g - mean) / std as doubles, so multiply becomes a pure dot product.
    void materialize_normalized_haploid();
    void materialize_normalized_diploid();

public:
    int start = 0;
    int end = 0;
    int num_samples = 0;
    int num_sites = 0;
    std::vector<int> positions;
    std::vector<ThreadingInstruction> instructions;

private:
    // Cached dense genotype matrix: genotype_matrix[site * num_samples + sample]
    std::vector<int> genotype_matrix;
    bool genotypes_materialized = false;

    // Cached diploid sum matrix: diploid_matrix[site * n_dip + i] = g[2i] + g[2i+1]
    std::vector<int> diploid_matrix;
    bool diploid_materialized = false;
    void materialize_diploid();

    // Cached standardized matrices for normalized multiply
    // standardized_hap[site * num_samples + i] = (g[i] - mu) / std
    std::vector<double> standardized_hap;
    bool standardized_hap_ready = false;
    // standardized_dip[site * n_dip + i] = ((g[2i]+g[2i+1]) - mu) / std
    std::vector<double> standardized_dip;
    bool standardized_dip_ready = false;

    // Tree-propagation multiply cache
    bool tree_ready = false;
    int tree_n_intervals = 0;
    std::vector<int> tree_ivl_start;   // first site index of each interval
    std::vector<int> tree_ivl_end;     // one past last site index
    std::vector<int> tree_ivl_seg;     // flat [sample * n_intervals + interval] -> segment
    std::vector<int> tree_ref_genome;

    // Precomputed mismatch correction signs: tree_mm_sign[offset + k] = +1 or -1
    // for non-self mismatches, 0 for self-ref (unused). Per-sample offsets in tree_mm_offset.
    std::vector<int8_t> tree_mm_sign;
    std::vector<int> tree_mm_offset;   // tree_mm_offset[i] = start index for sample i

    // Precomputed mismatch-to-interval mapping: tree_mm_ivl[offset + k] = interval index
    // containing mismatch k for sample i. Eliminates binary search during multiply.
    std::vector<int> tree_mm_ivl;

    // Precomputed carry_geno at end of each interval for each sample.
    // tree_carry[sample * n_intervals + interval] = genotype at last site of interval.
    std::vector<int8_t> tree_carry;

    // Per-sample segment-to-interval mapping for segment-level loop.
    // tree_seg_first_ivl[sample * max_segs + seg] = first interval index of segment.
    // Stored as flat vectors indexed by tree_seg_offset[sample] + seg.
    std::vector<int> tree_seg_first_ivl;  // first interval of each segment
    std::vector<int> tree_seg_offset;     // offset into tree_seg_first_ivl for each sample
};

#endif // THREADS_ARG_THREADING_INSTRUCTIONS_HPP
