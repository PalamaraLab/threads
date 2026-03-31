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

#include <functional>
#include <limits>
#include <set>
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

    // Add variants to an existing ThreadingInstructions.  Segments (starts,
    // targets, tmrcas) are unchanged — only positions and mismatches are updated.
    // new_positions: genomic coordinates of new variants (need not be sorted).
    // new_genotypes: row-major (n_new × num_samples), values 0 or 1.
    ThreadingInstructions add_variants(const std::vector<int>& new_positions,
                                       const std::vector<int>& new_genotypes,
                                       int n_new) const;

    // Simplify for multiply: strip tmrcas (set to 0) and merge consecutive
    // segments that share the same target. Returns a new ThreadingInstructions
    // with fewer segments. The genotype matrix is identical — only genealogical
    // information is discarded. Useful for multiply-only workloads.
    ThreadingInstructions simplify_for_multiply() const;

    // Aggressively coarsen by merging short segments into neighbors.
    // A segment covering fewer than min_sites variant sites is absorbed into
    // the previous segment (its target is replaced by the neighbor's target).
    // This is lossy — it changes some genotypes — but reduces tree_n_intervals.
    ThreadingInstructions coarsen(int min_sites) const;

    // Common operations
    std::vector<double> left_multiply(const std::vector<double>& x, bool diploid=false, bool normalize=false);
    std::vector<double> right_multiply(const std::vector<double>& x, bool diploid=false, bool normalize=false);

    // Tree-propagation multiply (tree shuttle).
    // Prepare: O(n*S*d) where S=segments/sample, d=tree depth.
    //   Builds interval decomposition, mismatch signs via lazy chain tracing,
    //   inverted mismatch index, interval change list, and ancestor chains.
    //   No genotype materialization; peak memory O(S*n + M + n*S*d).
    // Right multiply: O((n*ni + M)/T) per call, O(n*T) memory.
    //   Region-chunked OpenMP; each thread initializes via binary search.
    // Left multiply: O((m + M + S*d)/T) per call, O(n*T + m) memory.
    //   Sublinear via incremental O(n) weight vector (not O(n*ni) matrix).
    //   Region-chunked OpenMP; non-overlapping site ranges, no reduction.
    void prepare_tree_multiply();
    std::vector<double> right_multiply_tree(const std::vector<double>& x);
    std::vector<double> left_multiply_tree(const std::vector<double>& x);

    // Batch tree multiply: process k vectors in a single tree traversal.
    // Input/output are row-major flat arrays: X[row * k + col].
    // right_multiply_tree_batch: X is (num_sites, k), returns (num_samples, k)
    // left_multiply_tree_batch:  X is (num_samples, k), returns (num_sites, k)
    std::vector<double> right_multiply_tree_batch(const std::vector<double>& x_flat, int k);
    std::vector<double> left_multiply_tree_batch(const std::vector<double>& x_flat, int k);

    // RLE (run-length encoded) multiply: precomputes per-sample 1-runs,
    // then multiplies via prefix-sum lookups in O(m + total_1_runs).
    // Prepare is O(n*m) with O(m * tree_depth) peak memory (ref-counted cache).
    // Per-call cost is much lower than tree multiply when total_1_runs << n*ni.
    void prepare_rle_multiply();
    std::vector<double> right_multiply_rle(const std::vector<double>& x);
    std::vector<double> left_multiply_rle(const std::vector<double>& x);

    // Batch RLE multiply: process k vectors in one call.
    // Input/output are row-major flat arrays: X[row * k + col].
    // right: X is (num_sites, k), returns (num_samples, k)
    // left:  X is (num_samples, k), returns (num_sites, k)
    std::vector<double> right_multiply_rle_batch(const std::vector<double>& x_flat, int k);
    std::vector<double> left_multiply_rle_batch(const std::vector<double>& x_flat, int k);

    // DAG (mutation-set directed acyclic graph) multiply: preprocesses threading
    // instructions into a compact DAG of mutation sets, following the ARGMatMult
    // approach. Segments without mismatches are collapsed (pass-through), yielding
    // a DAG whose size grows sublinearly with n.
    //
    // The genomic region is split into num_chunks independent sub-DAGs (default:
    // OMP_NUM_THREADS). Each chunk covers a disjoint site range and can be
    // processed in parallel, giving embarrassingly-parallel speedup for the
    // propagation phase.
    //
    // Prepare: builds split points (reverse pass), then per-chunk DAG construction.
    // Right multiply: O(sum_c (|DAG_c| + connections_c) * k) — sublinear in n.
    // Left multiply:  O(sum_c (|DAG_c| + connections_c + m_c) * k) — sublinear in n.
    void prepare_dag_multiply(int num_chunks = 0);
    std::vector<double> right_multiply_dag(const std::vector<double>& x);
    std::vector<double> left_multiply_dag(const std::vector<double>& x);

    // Batch DAG multiply: process k vectors via outer-loop parallelism.
    // Each column is an independent multiply; OMP distributes columns across threads.
    // DAG structure is shared read-only in L3; per-thread data (partial, result) is private.
    // Input/output are row-major flat arrays: X[row * k + col].
    // right_multiply_dag_batch: X is (num_sites, k), returns (num_samples, k)
    // left_multiply_dag_batch:  X is (num_samples, k), returns (num_sites, k)
    std::vector<double> right_multiply_dag_batch(const std::vector<double>& x_flat, int k);
    std::vector<double> left_multiply_dag_batch(const std::vector<double>& x_flat, int k);

    // DAG statistics (valid after prepare_dag_multiply)
    int get_dag_num_sets() const;
    int get_dag_total_edges() const;
    int get_dag_total_connections() const;
    int get_dag_total_mutations() const;
    int get_dag_num_chunks() const { return dag_ready ? static_cast<int>(dag_chunks.size()) : 0; }

    // ARG traversals directly on threading instructions (no materialized ARG).
    //
    // visit_clades: enumerates all clades (internal nodes + leaves) in the
    // implicit threading tree, interval by interval.
    // callback(descendants, height, start_bp, end_bp)
    //   descendants: sorted vector of leaf sample indices under this clade
    //   height: coalescence height of the clade (0 for leaves)
    //   start_bp, end_bp: physical position range
    void visit_clades(std::function<void(const std::vector<int>&, double, int, int)> callback) const;

    // visit_branches: enumerates all branches (edges between nodes).
    // callback(descendants, child_height, parent_height, start_bp, end_bp)
    //   descendants: sorted vector of leaf sample indices below this branch
    //   child_height: height of the child node (0 for leaf branches)
    //   parent_height: height of the parent node
    void visit_branches(std::function<void(const std::vector<int>&, double, double, int, int)> callback) const;

    // ARG summary statistics (optimized, no full descendant tracking).
    //
    // total_volume: sum of (parent_height - child_height) * (end_bp - start_bp).
    double total_volume() const;
    //
    // allele_frequency_spectrum: afs[k] = total branch volume for branches
    //   with exactly k descendants, for k in 0..num_samples.
    std::vector<double> allele_frequency_spectrum() const;
    //
    // association_diploid: for each clade, compute chi-square association
    //   statistic against a diploid phenotype vector. Returns (positions,
    //   chi2_stats, p_values, n_descendants) for all clades with p < threshold.
    //   phenotypes has length num_samples/2.
    struct AssociationResult {
        std::vector<int> position;         // midpoint bp
        std::vector<double> chi2;          // chi-square statistic
        std::vector<double> pvalue;        // p-value
        std::vector<int> n_descendants;    // clade size (haploid)
    };
    AssociationResult association_diploid(const std::vector<double>& phenotypes,
                                          double p_threshold = 5e-8) const;
    //
    // generate_mutations: Poisson-sample mutations on branches proportional to
    //   branch volume. Returns (positions_bp, genotype_flat) where genotype_flat
    //   is row-major (n_mutations × num_samples) with 0/1 values.
    struct MutationResult {
        std::vector<int> positions;
        std::vector<int> genotypes;  // flat row-major, n_mutations * num_samples
        int n_mutations;
    };
    MutationResult generate_mutations(double mutation_rate, int seed = 42) const;

    // Count heterozygous diploid individuals per site without materializing
    // the full genotype matrix.  Returns vector of length num_sites where
    // het[s] = number of diploid pairs (2j, 2j+1) with different alleles.
    // O(n*m) time, O(n) memory (no n*m matrix allocated).
    std::vector<int> het_per_site() const;

    // Count heterozygous sites per diploid individual.  Returns vector of
    // length num_samples/2 where het[j] = number of sites where sample 2j
    // and sample 2j+1 have different alleles.  O(n*m) time, O(n) memory.
    std::vector<int> het_per_individual() const;

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

    // Returns the dense genotype matrix (materializes if needed).
    const std::vector<int>& get_genotype_matrix() { materialize_genotypes(); return genotype_matrix; }

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
    std::vector<int> tree_ref_genome;

    // Precomputed mismatch correction signs: tree_mm_sign[offset + k] = +1 or -1
    // for non-self mismatches, 0 for self-ref (unused). Per-sample offsets in tree_mm_offset.
    std::vector<int8_t> tree_mm_sign;
    std::vector<int> tree_mm_offset;   // tree_mm_offset[i] = start index for sample i

    // Site-to-interval lookup: site_to_ivl[s] = interval index containing site s.
    // O(m) memory, replaces O(M_total) tree_mm_ivl array.
    std::vector<int> site_to_ivl;

    // Per-sample segment-to-interval mapping for segment-level loop.
    // tree_seg_first_ivl[offset + k] = first global interval of k-th segment
    // tree_seg_id[offset + k] = segment index in instructions[i].targets
    // Indexed by tree_seg_offset[sample] + k.
    std::vector<int> tree_seg_first_ivl;
    std::vector<int> tree_seg_id;         // segment index for target lookup
    std::vector<int> tree_seg_offset;

    // Lazy genotype query: trace threading chain to sample 0.
    // O(d * log S) per call, d = tree depth.
    int genotype_at_site(int sample, int site) const;

    // Target lookup: which sample does `sample` target at a given interval?
    // O(log S_i) via binary search on tree_seg_first_ivl.
    int target_at_interval(int sample, int interval) const;

    // Inverted mismatch index: non-self-ref mismatches grouped by interval.
    // inv_mm_offset[j]..inv_mm_offset[j+1] indexes into the flat arrays.
    std::vector<int> inv_mm_sample;
    std::vector<int> inv_mm_site;
    std::vector<int8_t> inv_mm_sign;
    std::vector<int> inv_mm_offset;    // ni+1 entries

    // Interval change list: samples whose target changes at each boundary.
    // Sorted by sample descending within each boundary group.
    std::vector<int> ivl_change_sample;
    std::vector<int> ivl_change_old_tgt;
    std::vector<int> ivl_change_new_tgt;
    std::vector<int> ivl_change_offset;  // ni+1 entries

    // Self-ref samples per interval (for carry-run handling in left multiply).
    std::vector<int> self_ref_sample;
    std::vector<int> self_ref_offset;    // ni+1 entries

    // Precomputed ancestor chains for left multiply W updates.
    // chain_data[chain_offset[c]..chain_offset[c+1]) = ancestor list
    // for change c (detach uses old tree, attach uses new tree).
    std::vector<int> detach_chain_data;
    std::vector<int> detach_chain_offset;  // total_changes + 1
    std::vector<int> attach_chain_data;
    std::vector<int> attach_chain_offset;  // total_changes + 1


    // RLE multiply cache: per-sample 1-runs stored as flat (start, end) pairs.
    // rle_run_start[rle_offset[i] .. rle_offset[i+1]) = start site indices of 1-runs
    // rle_run_end  [rle_offset[i] .. rle_offset[i+1]) = end site indices of 1-runs
    bool rle_ready = false;
    std::vector<int> rle_offset;      // n+1 entries
    std::vector<int> rle_run_start;   // total_runs entries
    std::vector<int> rle_run_end;     // total_runs entries

    // DAG multiply cache: region-chunked mutation-set DAGs.
    struct DagChunk {
        int num_sets = 0;
        int site_lo = 0, site_hi = 0;  // site index range [lo, hi)
        std::vector<int> d2a, d2a_off;
        std::vector<int> a2d, a2d_off;
        std::vector<int> mut_site;
        std::vector<int8_t> mut_sign;
        std::vector<int> mut_off;
        std::vector<int> indiv_set, indiv_off;
        std::vector<int> set_indiv, set_indiv_off;
        std::vector<int> ref_ones;
    };
    bool dag_ready = false;
    std::vector<DagChunk> dag_chunks;

    DagChunk build_dag_chunk(const std::vector<std::set<int>>& sample_splits,
                             int chunk_site_lo, int chunk_site_hi);

    // Non-OMP single-column multiply helpers (called within OMP parallel regions)
    void right_multiply_dag_single(const double* x, double* out) const;
    void left_multiply_dag_single(const double* y, double* out) const;
};

#endif // THREADS_ARG_THREADING_INSTRUCTIONS_HPP
