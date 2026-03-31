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

#include "ThreadingInstructions.hpp"
#include "GenotypeIterator.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>
#include <random>
#include <iostream>
#include <limits>
#include <deque>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


// perf: Note that args are passed by value rather than ref as they run faster
// when used with std::move below.
ThreadingInstruction::ThreadingInstruction(std::vector<int> _starts,
                                           std::vector<double> _tmrcas,
                                           std::vector<int> _targets,
                                           std::vector<int> _mismatches) :
    starts(_starts), tmrcas(_tmrcas), targets(_targets), mismatches(_mismatches) {
    num_segments = starts.size();
    num_mismatches = mismatches.size();
    if (tmrcas.size() != num_segments || targets.size() != num_segments) {
      throw std::runtime_error("Mismatching lengths of threading instruction input");
    }
}

ThreadingInstructionIterator::ThreadingInstructionIterator(
    const ThreadingInstruction& _instruction,
    const std::vector<int>& _positions
) :
    instruction(_instruction), positions(_positions) {
    if (instruction.num_segments == 0) {
        current_tmrca = 0;
        current_target = -1;
    } else {
        current_tmrca = instruction.tmrcas.front();
        current_target = instruction.targets.front();
    }
    if (instruction.num_segments <= 1) {
        next_segment_start = std::numeric_limits<int>::max();
    } else {
        next_segment_start = instruction.starts.at(1);
    }

    if (instruction.num_mismatches == 0) {
        next_mismatch = std::numeric_limits<int>::max();
        is_mismatch = false;
    } else {
        next_mismatch = positions.at(instruction.mismatches.front());
        is_mismatch = instruction.mismatches.front() == 0;
    }
}

void ThreadingInstructionIterator::increment_site(int position) {
    // Search forward and update current target/tmrca for this position
    while (next_segment_start <= position && next_segment_start != std::numeric_limits<int>::max()) {
        current_segment++;
        current_target = instruction.targets.at(current_segment);
        current_tmrca = instruction.tmrcas.at(current_segment);
        if (current_segment < instruction.num_segments - 1) {
            next_segment_start = instruction.starts.at(current_segment + 1);
        } else {
            next_segment_start = std::numeric_limits<int>::max();
        }
    }

    // Search forward and update current mismatch status
    while (next_mismatch < position ) {
        current_mismatch++;
        if (current_mismatch <= instruction.num_mismatches - 1) {
            next_mismatch = positions.at(instruction.mismatches.at(current_mismatch));
        } else {
            next_mismatch = std::numeric_limits<int>::max();
        }
    }
    is_mismatch = position == next_mismatch;

    sites_processed++;
}

ThreadingInstructions::ThreadingInstructions(std::vector<ThreadingInstruction> _instructions, std::vector<int> _positions) :
    instructions(_instructions), positions(_positions) {
    start = positions.front();
    end = positions.back() + 1;
    num_samples = instructions.size();
    num_sites = positions.size();
}

ThreadingInstructions::ThreadingInstructions(const std::vector<ViterbiPath>& paths, const int _start, const int _end, const std::vector<int>& all_positions) :
    start(_start), end(_end) {
    num_samples = paths.size();

    int pos_idx_offset = -1;
    int i = 0;
    for (auto pos : all_positions) {
        if (start <= pos && pos < end) {
            positions.push_back(pos);
            pos_idx_offset = pos_idx_offset < 0 ? i : pos_idx_offset;
        }
        if (pos > end) {
            break;
        }
        i++;
    }
    pos_idx_offset = pos_idx_offset < 0 ? 0 : pos_idx_offset;
    num_sites = positions.size();

    std::vector<int> starts;
    std::vector<int> targets;
    std::vector<double> tmrcas;
    instructions.reserve(num_samples);
    for (auto path : paths) {
        std::vector<int> mismatches;
        path.map_positions(all_positions);
        std::tie(starts, targets, tmrcas) = path.dump_data_in_range(start, end);
        for (auto het_idx : path.het_sites) {
            int het_pos = all_positions.at(het_idx);
            if ((start <= het_pos) && (het_pos < end)) {
                mismatches.push_back(het_idx - pos_idx_offset);
            }
        }
        instructions.emplace_back(starts, tmrcas, targets, mismatches);
    }
}

ThreadingInstructions::ThreadingInstructions(const std::vector<std::vector<int>>& starts,
                        const std::vector<std::vector<double>>& tmrcas,
                        const std::vector<std::vector<int>>& targets,
                        const std::vector<std::vector<int>>& mismatches,
                        const std::vector<int>& _positions, int _start, int _end) :
    positions(_positions), start(_start), end(_end) {
    num_samples = starts.size();
    num_sites = positions.size();
    instructions.reserve(num_samples);
    for (int i = 0; i < num_samples; i++) {
        instructions.emplace_back(starts.at(i), tmrcas.at(i), targets.at(i), mismatches.at(i));
    }
}

std::vector<std::vector<int>> ThreadingInstructions::all_starts() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.starts);
    }
    return out;
}

std::vector<std::vector<double>> ThreadingInstructions::all_tmrcas() {
    std::vector<std::vector<double>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.tmrcas);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_targets() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.targets);
    }
    return out;
}

std::vector<std::vector<int>> ThreadingInstructions::all_mismatches() {
    std::vector<std::vector<int>> out;
    out.reserve(num_samples);
    for (auto& instruction : instructions) {
        out.push_back(instruction.mismatches);
    }
    return out;
}

ThreadingInstructions ThreadingInstructions::sub_range(const int range_start, const int range_end) const {
    // Validate range and get corresponding start end index in positions
    if (range_start > range_end) {
        std::ostringstream oss;
        oss << "range_start " << range_start << " cannot be greater than range_end " << range_end;
        throw std::runtime_error(oss.str());
    }

    const auto range_start_it = std::lower_bound(positions.begin(), positions.end(), range_start);
    if ((range_start < 0) || (range_start_it == positions.end())) {
        std::ostringstream oss;
        oss << "range_start " << range_start << " is not a valid position";
        throw std::runtime_error(oss.str());
    }
    const int range_start_idx = range_start_it - positions.begin();

    const auto range_end_it = std::lower_bound(positions.begin(), positions.end(), range_end);
    if ((range_end < 0) || (range_end_it == positions.end())) {
        std::ostringstream oss;
        oss << "range_end " << range_end << " is not a valid position";
        throw std::runtime_error(oss.str());
    }
    const int range_end_idx = range_end_it - positions.begin();

    // Create range positions from indexed subset
    std::vector<int> range_positions{
        positions.begin() + range_start_idx,
        positions.begin() + range_end_idx + 1
    };

    std::vector<ThreadingInstruction> range_instructions;
    for (const auto& instruction : instructions) {
        std::vector<int> sub_starts;
        std::vector<double> sub_tmrcas;
        std::vector<int> sub_targets;
        std::vector<int> sub_mismatches;

        // Create sub-vectors based on range
        const size_t num_starts = instruction.starts.size();
        for (size_t pos_idx = 0; pos_idx < num_starts; ++pos_idx) {
            const int seg_start = instruction.starts[pos_idx];
            const int seg_end = (pos_idx == num_starts - 1) ?
                std::numeric_limits<int>::max() :
                instruction.starts[pos_idx + 1];

            if ((range_start <= seg_end) && (range_end >= seg_start)) {
                sub_starts.push_back(seg_start);
                sub_tmrcas.push_back(instruction.tmrcas[pos_idx]);
                sub_targets.push_back(instruction.targets[pos_idx]);
            }
        }

        // Recompute mismatches based on start indexes within range
        for (const int mismatch_idx : instruction.mismatches) {
            const int pos = positions[mismatch_idx];
            if ((pos >= range_start) && (pos <= range_end)) {
                sub_mismatches.push_back(mismatch_idx - range_start_idx);
            }
        }

        if (!sub_starts.empty() || !sub_mismatches.empty()) {
            range_instructions.emplace_back(
                std::move(sub_starts),
                std::move(sub_tmrcas),
                std::move(sub_targets),
                std::move(sub_mismatches)
            );
        }
    }

    return ThreadingInstructions{
        std::move(range_instructions),
        std::move(range_positions)
    };
}

ThreadingInstructions ThreadingInstructions::add_variants(
        const std::vector<int>& new_positions_unsorted,
        const std::vector<int>& new_genotypes_flat,
        int n_new) const {
    const int n_old = num_sites;
    const int n_samp = num_samples;

    // Sort new positions by value, keeping a permutation for genotype reorder
    std::vector<int> order(n_new);
    std::iota(order.begin(), order.end(), 0);
    std::vector<int> new_pos(n_new);
    for (int i = 0; i < n_new; i++) new_pos[i] = new_positions_unsorted[i];
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return new_pos[a] < new_pos[b]; });
    {
        std::vector<int> sorted_pos(n_new);
        for (int i = 0; i < n_new; i++) sorted_pos[i] = new_pos[order[i]];
        new_pos = std::move(sorted_pos);
    }

    // Merge old and new positions; record index mappings
    std::vector<int> merged_positions;
    merged_positions.reserve(n_old + n_new);
    std::vector<int> old_to_merged(n_old);
    std::vector<int> new_to_merged(n_new);

    int oi = 0, ni = 0;
    while (oi < n_old && ni < n_new) {
        if (positions[oi] <= new_pos[ni]) {
            old_to_merged[oi] = static_cast<int>(merged_positions.size());
            merged_positions.push_back(positions[oi]);
            oi++;
        } else {
            new_to_merged[ni] = static_cast<int>(merged_positions.size());
            merged_positions.push_back(new_pos[ni]);
            ni++;
        }
    }
    while (oi < n_old) {
        old_to_merged[oi] = static_cast<int>(merged_positions.size());
        merged_positions.push_back(positions[oi]);
        oi++;
    }
    while (ni < n_new) {
        new_to_merged[ni] = static_cast<int>(merged_positions.size());
        merged_positions.push_back(new_pos[ni]);
        ni++;
    }

    // Site-major processing: for each new site, compute all samples' genotypes.
    // Memory: O(n_samp) for the current site's genotypes + O(n_samp) for segment
    // cursors. Total working memory: O(n_samp) regardless of n_new.
    //
    // Each sample maintains a segment cursor that advances monotonically since
    // new sites are sorted. Cost per site per sample: O(1) amortized.
    std::vector<std::vector<int>> per_sample_new_mm(n_samp);
    std::vector<int8_t> site_geno(n_samp);     // genotypes at current site
    std::vector<int> seg_cursor(n_samp, 0);     // per-sample segment cursor

    for (int j = 0; j < n_new; j++) {
        const int pos = new_pos[j];
        const size_t src_row = static_cast<size_t>(order[j]) * n_samp;

        for (int s = 0; s < n_samp; s++) {
            int8_t my_geno = static_cast<int8_t>(new_genotypes_flat[src_row + s]);

            if (s == 0) {
                site_geno[0] = my_geno;
                if (my_geno == 1) {
                    per_sample_new_mm[0].push_back(new_to_merged[j]);
                }
            } else {
                const auto& seg_starts = instructions[s].starts;
                const auto& seg_targets = instructions[s].targets;
                const int n_segs = static_cast<int>(seg_starts.size());

                // Advance cursor (amortized O(1) across all sites)
                int& cur = seg_cursor[s];
                while (cur + 1 < n_segs && seg_starts[cur + 1] <= pos) {
                    cur++;
                }

                int target_geno = site_geno[seg_targets[cur]];
                site_geno[s] = my_geno;

                if (my_geno != target_geno) {
                    per_sample_new_mm[s].push_back(new_to_merged[j]);
                }
            }
        }
    }

    // Build final instructions: merge old (remapped) and new mismatch lists
    std::vector<ThreadingInstruction> new_instructions;
    new_instructions.reserve(n_samp);

    for (int s = 0; s < n_samp; s++) {
        const auto& inst = instructions[s];

        // Remap old mismatches to merged indices (preserves sorted order)
        std::vector<int> old_remapped;
        old_remapped.reserve(inst.mismatches.size());
        for (int m : inst.mismatches) {
            old_remapped.push_back(old_to_merged[m]);
        }

        // Merge two sorted sequences
        const auto& nm = per_sample_new_mm[s];
        std::vector<int> merged_mm(old_remapped.size() + nm.size());
        std::merge(old_remapped.begin(), old_remapped.end(),
                   nm.begin(), nm.end(),
                   merged_mm.begin());

        new_instructions.emplace_back(
            std::vector<int>(inst.starts),
            std::vector<double>(inst.tmrcas),
            std::vector<int>(inst.targets),
            std::move(merged_mm)
        );
    }

    ThreadingInstructions result(std::move(new_instructions), std::move(merged_positions));
    result.start = start;
    result.end = end;
    return result;
}

ThreadingInstructions ThreadingInstructions::simplify_for_multiply() const {
    std::vector<ThreadingInstruction> simplified;
    simplified.reserve(num_samples);

    for (const auto& inst : instructions) {
        std::vector<int> new_starts;
        std::vector<double> new_tmrcas;
        std::vector<int> new_targets;

        if (inst.num_segments > 0) {
            new_starts.push_back(inst.starts[0]);
            new_tmrcas.push_back(0.0);
            new_targets.push_back(inst.targets[0]);

            for (size_t k = 1; k < inst.num_segments; k++) {
                if (inst.targets[k] != new_targets.back()) {
                    new_starts.push_back(inst.starts[k]);
                    new_tmrcas.push_back(0.0);
                    new_targets.push_back(inst.targets[k]);
                }
                // else: same target, merge by skipping this segment boundary
            }
        }

        // Mismatches are unchanged — they index into positions[], not segments
        simplified.emplace_back(
            std::move(new_starts),
            std::move(new_tmrcas),
            std::move(new_targets),
            std::vector<int>(inst.mismatches)
        );
    }

    ThreadingInstructions result(std::move(simplified), std::vector<int>(positions));
    result.start = start;
    result.end = end;
    return result;
}

ThreadingInstructions ThreadingInstructions::coarsen(int min_sites) const {
    // For each sample, merge segments that cover fewer than min_sites variant
    // sites into their predecessor. Mismatches are recomputed for altered
    // segments by comparing the sample's true genotype against the new target.

    // Precompute reference genome (sample 0 mismatches)
    std::vector<int8_t> ref_genome(num_sites, 0);
    for (int idx : instructions[0].mismatches) {
        if (idx >= 0 && idx < num_sites) ref_genome[idx] = 1;
    }

    // Build full genotype matrix once (needed to recompute mismatches after
    // target changes). Use the same tree-propagation as GenotypeIterator.
    // geno[sample * num_sites + site]
    std::vector<int8_t> geno(static_cast<size_t>(num_samples) * num_sites);

    // Sample 0
    for (int s = 0; s < num_sites; s++) {
        geno[s] = ref_genome[s];
    }

    // Samples 1..n-1
    for (int i = 1; i < num_samples; i++) {
        const auto& inst = instructions[i];
        int seg = 0;
        int mm_idx = 0;
        const int n_segs = static_cast<int>(inst.starts.size());
        const int n_mm = static_cast<int>(inst.mismatches.size());
        const size_t row = static_cast<size_t>(i) * num_sites;

        for (int s = 0; s < num_sites; s++) {
            // Advance segment
            while (seg + 1 < n_segs && inst.starts[seg + 1] <= positions[s]) {
                seg++;
            }
            int target = inst.targets[seg];
            int tgt_g = geno[static_cast<size_t>(target) * num_sites + s];
            bool is_mm = (mm_idx < n_mm && inst.mismatches[mm_idx] == s);
            if (is_mm) mm_idx++;
            geno[row + s] = static_cast<int8_t>(is_mm ? (1 - tgt_g) : tgt_g);
        }
    }

    // Now coarsen each sample's segments
    std::vector<ThreadingInstruction> new_instructions;
    new_instructions.reserve(num_samples);

    for (int i = 0; i < num_samples; i++) {
        const auto& inst = instructions[i];
        const int n_segs = static_cast<int>(inst.starts.size());

        if (i == 0 || n_segs <= 1) {
            // Sample 0 or single-segment: keep as-is
            new_instructions.emplace_back(
                std::vector<int>(inst.starts),
                std::vector<double>(inst.tmrcas),
                std::vector<int>(inst.targets),
                std::vector<int>(inst.mismatches));
            continue;
        }

        // Compute site count per segment
        std::vector<int> seg_site_count(n_segs, 0);
        {
            int seg = 0;
            for (int s = 0; s < num_sites; s++) {
                while (seg + 1 < n_segs && inst.starts[seg + 1] <= positions[s]) {
                    seg++;
                }
                seg_site_count[seg]++;
            }
        }

        // Merge short segments into predecessor
        std::vector<int> new_starts, new_targets;
        std::vector<double> new_tmrcas;
        new_starts.push_back(inst.starts[0]);
        new_targets.push_back(inst.targets[0]);
        new_tmrcas.push_back(0.0);

        for (int k = 1; k < n_segs; k++) {
            if (seg_site_count[k] < min_sites) {
                // Absorb into previous segment (skip this boundary)
                continue;
            }
            new_starts.push_back(inst.starts[k]);
            new_targets.push_back(inst.targets[k]);
            new_tmrcas.push_back(0.0);
        }

        // Recompute mismatches with the new target assignments
        std::vector<int> new_mm;
        const size_t row = static_cast<size_t>(i) * num_sites;
        int seg = 0;
        const int new_n_segs = static_cast<int>(new_starts.size());

        for (int s = 0; s < num_sites; s++) {
            while (seg + 1 < new_n_segs && new_starts[seg + 1] <= positions[s]) {
                seg++;
            }
            int target = new_targets[seg];
            int tgt_g = geno[static_cast<size_t>(target) * num_sites + s];
            int my_g = geno[row + s];
            if (my_g != tgt_g) {
                new_mm.push_back(s);
            }
        }

        new_instructions.emplace_back(
            std::move(new_starts), std::move(new_tmrcas),
            std::move(new_targets), std::move(new_mm));
    }

    ThreadingInstructions result(std::move(new_instructions), std::vector<int>(positions));
    result.start = start;
    result.end = end;
    return result;
}

// ── RLE multiply ────────────────────────────────────────────────────────────
//
// RLE multiply is equivalent to a sparse matrix representation: each sample's
// genotype is stored as runs of 1s (start, end pairs). Multiply accumulates
// contributions from these runs via prefix sums.
//
// Complexity: O(m + total_1_runs) per multiply, where total_1_runs is the sum
// of 1-run counts across all samples. In practice this is O(n * m * p̄) where
// p̄ is the average allele frequency — the same dominant term as tree shuttle's
// mismatch cost. However, tree shuttle replaces the O(m) prefix-sum sweep with
// O(n * n_intervals) tree propagation, which is typically 10-100x smaller than
// m, so tree shuttle is preferred for large m.

void ThreadingInstructions::prepare_rle_multiply() {
    if (rle_ready) return;

    const int n = num_samples;
    const int m = num_sites;

    rle_offset.resize(n + 1);
    std::vector<int> run_starts, run_ends;

    // Reference counting for genotype cache (same as tree prepare)
    std::vector<int> ref_count(n, 0);
    for (int i = 1; i < n; i++) {
        std::set<int> seen;
        for (int t : instructions[i].targets) {
            if (t != i && seen.insert(t).second) ref_count[t]++;
        }
    }

    std::unordered_map<int, std::vector<uint8_t>> geno_cache;
    std::vector<uint8_t> geno_buf(m, 0);

    for (int i = 0; i < n; i++) {
        const auto& inst_i = instructions[i];
        const int n_segs = static_cast<int>(inst_i.num_segments);
        const int* mm_data = inst_i.mismatches.data();
        const int n_mm = static_cast<int>(inst_i.num_mismatches);

        if (i == 0) {
            // Sample 0: genotype from reference mismatches
            std::fill(geno_buf.begin(), geno_buf.end(), 0);
            for (int k = 0; k < n_mm; k++) {
                if (mm_data[k] >= 0 && mm_data[k] < m) geno_buf[mm_data[k]] = 1;
            }
        } else {
            // Build genotype from segments
            for (int seg = 0; seg < n_segs; seg++) {
                // Convert segment start position to site index
                int site_start = 0;
                if (seg > 0) {
                    auto it = std::lower_bound(positions.begin(), positions.end(),
                                               inst_i.starts[seg]);
                    site_start = static_cast<int>(it - positions.begin());
                }
                int site_end = m;
                if (seg + 1 < n_segs) {
                    auto it = std::lower_bound(positions.begin(), positions.end(),
                                               inst_i.starts[seg + 1]);
                    site_end = static_cast<int>(it - positions.begin());
                }

                const int target = inst_i.targets[seg];

                if (target != i) {
                    const auto& tgt = geno_cache.at(target);
                    std::memcpy(&geno_buf[site_start], &tgt[site_start],
                                site_end - site_start);
                    // Flip at mismatches
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, site_start);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, site_end);
                    for (const int* it = lo; it != hi; ++it) {
                        geno_buf[*it] ^= 1;
                    }
                } else {
                    // Self-ref: carry forward, flip at mismatches
                    int carry = (site_start > 0) ? geno_buf[site_start - 1] : 0;
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, site_start);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, site_end);
                    const int* it = lo;
                    for (int s = site_start; s < site_end; s++) {
                        if (it != hi && *it == s) {
                            carry ^= 1;
                            ++it;
                        }
                        geno_buf[s] = static_cast<uint8_t>(carry);
                    }
                }
            }
        }

        // Extract 1-runs
        rle_offset[i] = static_cast<int>(run_starts.size());
        bool in_run = false;
        for (int s = 0; s < m; s++) {
            if (geno_buf[s] && !in_run) {
                run_starts.push_back(s);
                in_run = true;
            } else if (!geno_buf[s] && in_run) {
                run_ends.push_back(s);
                in_run = false;
            }
        }
        if (in_run) run_ends.push_back(m);

        // Cache if still referenced
        if (ref_count[i] > 0) {
            geno_cache[i] = geno_buf;
        }
        // Decrement refs
        if (i > 0) {
            std::set<int> seen;
            for (int t : inst_i.targets) {
                if (t != i && seen.insert(t).second) {
                    if (--ref_count[t] == 0) geno_cache.erase(t);
                }
            }
        }
    }
    rle_offset[n] = static_cast<int>(run_starts.size());

    rle_run_start = std::move(run_starts);
    rle_run_end = std::move(run_ends);
    rle_ready = true;
}

std::vector<double> ThreadingInstructions::right_multiply_rle(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_sites) {
        throw std::runtime_error("Input vector must have length num_sites.");
    }
    prepare_rle_multiply();

    const int n = num_samples;
    const int m = num_sites;

    // Prefix sums of x
    std::vector<double> prefix(m + 1, 0.0);
    for (int s = 0; s < m; s++) prefix[s + 1] = prefix[s] + x[s];

    std::vector<double> out(n);
    const double* pfx = prefix.data();
    const int* rs_data = rle_run_start.data();
    const int* re_data = rle_run_end.data();

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        const int rs = rle_offset[i];
        const int re = rle_offset[i + 1];
        double sum = 0.0;
        for (int r = rs; r < re; r++) {
            sum += pfx[re_data[r]] - pfx[rs_data[r]];
        }
        out[i] = sum;
    }
    return out;
}

std::vector<double> ThreadingInstructions::left_multiply_rle(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_samples) {
        throw std::runtime_error("Input vector must have length num_samples.");
    }
    prepare_rle_multiply();

    const int n = num_samples;
    const int m = num_sites;

#ifdef _OPENMP
    const int nt = omp_get_max_threads();
    std::vector<double> all_diff(static_cast<size_t>(m + 1) * nt, 0.0);

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double* my_diff = &all_diff[static_cast<size_t>(tid) * (m + 1)];

        #pragma omp for schedule(static)
        for (int i = 0; i < n; i++) {
            const int rs = rle_offset[i];
            const int re = rle_offset[i + 1];
            const double xi = x[i];
            for (int r = rs; r < re; r++) {
                my_diff[rle_run_start[r]] += xi;
                my_diff[rle_run_end[r]] -= xi;
            }
        }
    }

    // Reduce thread-local diffs and prefix sum
    std::vector<double> out(m);
    double running = 0.0;
    for (int s = 0; s < m; s++) {
        double d = 0.0;
        for (int t = 0; t < nt; t++) {
            d += all_diff[static_cast<size_t>(t) * (m + 1) + s];
        }
        running += d;
        out[s] = running;
    }
#else
    // Single-threaded: diff array approach
    std::vector<double> diff(m + 1, 0.0);
    for (int i = 0; i < n; i++) {
        const int rs = rle_offset[i];
        const int re = rle_offset[i + 1];
        const double xi = x[i];
        for (int r = rs; r < re; r++) {
            diff[rle_run_start[r]] += xi;
            diff[rle_run_end[r]] -= xi;
        }
    }

    std::vector<double> out(m);
    double running = 0.0;
    for (int s = 0; s < m; s++) {
        running += diff[s];
        out[s] = running;
    }
#endif
    return out;
}

std::vector<double> ThreadingInstructions::right_multiply_rle_batch(
        const std::vector<double>& x_flat, int k) {
    if (k == 1) return right_multiply_rle(x_flat);
    if (static_cast<int>(x_flat.size()) != num_sites * k) {
        throw std::runtime_error("Input must have length num_sites * k.");
    }
    prepare_rle_multiply();

    const int n = num_samples;
    const int m = num_sites;

    // Prefix sums for each of k vectors: prefix[(s+1)*k + j] = sum x[0..s, j]
    std::vector<double> prefix(static_cast<size_t>(m + 1) * k, 0.0);
    for (int s = 0; s < m; s++) {
        const double* xrow = &x_flat[static_cast<size_t>(s) * k];
        const double* prev = &prefix[static_cast<size_t>(s) * k];
        double* cur = &prefix[static_cast<size_t>(s + 1) * k];
        for (int j = 0; j < k; j++) {
            cur[j] = prev[j] + xrow[j];
        }
    }

    std::vector<double> out(static_cast<size_t>(n) * k, 0.0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        const int rs = rle_offset[i];
        const int re = rle_offset[i + 1];
        double* out_i = &out[static_cast<size_t>(i) * k];
        for (int r = rs; r < re; r++) {
            const double* p_end = &prefix[static_cast<size_t>(rle_run_end[r]) * k];
            const double* p_start = &prefix[static_cast<size_t>(rle_run_start[r]) * k];
            for (int j = 0; j < k; j++) {
                out_i[j] += p_end[j] - p_start[j];
            }
        }
    }
    return out;
}

std::vector<double> ThreadingInstructions::left_multiply_rle_batch(
        const std::vector<double>& x_flat, int k) {
    if (k == 1) return left_multiply_rle(x_flat);
    if (static_cast<int>(x_flat.size()) != num_samples * k) {
        throw std::runtime_error("Input must have length num_samples * k.");
    }
    prepare_rle_multiply();

    const int n = num_samples;
    const int m = num_sites;

    const size_t diff_len = static_cast<size_t>(m + 1) * k;

#ifdef _OPENMP
    const int nt = omp_get_max_threads();
    std::vector<double> all_diff(diff_len * nt, 0.0);

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double* my_diff = &all_diff[static_cast<size_t>(tid) * diff_len];

        #pragma omp for schedule(static)
        for (int i = 0; i < n; i++) {
            const int rs = rle_offset[i];
            const int re = rle_offset[i + 1];
            const double* x_i = &x_flat[static_cast<size_t>(i) * k];
            for (int r = rs; r < re; r++) {
                double* d_start = &my_diff[static_cast<size_t>(rle_run_start[r]) * k];
                double* d_end = &my_diff[static_cast<size_t>(rle_run_end[r]) * k];
                for (int j = 0; j < k; j++) {
                    d_start[j] += x_i[j];
                    d_end[j] -= x_i[j];
                }
            }
        }
    }

    // Reduce + prefix sum
    std::vector<double> out(static_cast<size_t>(m) * k);
    std::vector<double> running(k, 0.0);
    for (int s = 0; s < m; s++) {
        double* o = &out[static_cast<size_t>(s) * k];
        for (int j = 0; j < k; j++) {
            double d = 0.0;
            for (int t = 0; t < nt; t++) {
                d += all_diff[static_cast<size_t>(t) * diff_len + static_cast<size_t>(s) * k + j];
            }
            running[j] += d;
            o[j] = running[j];
        }
    }
#else
    std::vector<double> diff(diff_len, 0.0);
    for (int i = 0; i < n; i++) {
        const int rs = rle_offset[i];
        const int re = rle_offset[i + 1];
        const double* x_i = &x_flat[static_cast<size_t>(i) * k];
        for (int r = rs; r < re; r++) {
            double* d_start = &diff[static_cast<size_t>(rle_run_start[r]) * k];
            double* d_end = &diff[static_cast<size_t>(rle_run_end[r]) * k];
            for (int j = 0; j < k; j++) {
                d_start[j] += x_i[j];
                d_end[j] -= x_i[j];
            }
        }
    }

    std::vector<double> out(static_cast<size_t>(m) * k);
    std::vector<double> running(k, 0.0);
    for (int s = 0; s < m; s++) {
        const double* d = &diff[static_cast<size_t>(s) * k];
        double* o = &out[static_cast<size_t>(s) * k];
        for (int j = 0; j < k; j++) {
            running[j] += d[j];
            o[j] = running[j];
        }
    }
#endif
    return out;
}

std::vector<int> ThreadingInstructions::het_per_site() const {
    const int n = num_samples;
    const int n_dip = n / 2;
    std::vector<int> het(num_sites, 0);

    GenotypeIterator gi(*this);
    int s = 0;
    while (gi.has_next_genotype()) {
        const std::vector<int>& g = gi.next_genotype();
        int h = 0;
        for (int j = 0; j < n_dip; j++) {
            h += (g[2 * j] != g[2 * j + 1]);
        }
        het[s++] = h;
    }
    return het;
}

std::vector<int> ThreadingInstructions::het_per_individual() const {
    const int n = num_samples;
    const int n_dip = n / 2;
    std::vector<int> het(n_dip, 0);

    GenotypeIterator gi(*this);
    while (gi.has_next_genotype()) {
        const std::vector<int>& g = gi.next_genotype();
        for (int j = 0; j < n_dip; j++) {
            het[j] += (g[2 * j] != g[2 * j + 1]);
        }
    }
    return het;
}

void ThreadingInstructions::materialize_genotypes() {
    if (genotypes_materialized) return;
    genotype_matrix.resize(static_cast<size_t>(num_sites) * num_samples);
    GenotypeIterator gi(*this);
    int site = 0;
    while (gi.has_next_genotype()) {
        const std::vector<int>& g = gi.next_genotype();
        std::copy(g.begin(), g.end(), genotype_matrix.begin() + static_cast<size_t>(site) * num_samples);
        site++;
    }
    genotypes_materialized = true;
}

void ThreadingInstructions::materialize_diploid() {
    if (diploid_materialized) return;
    materialize_genotypes();
    const int n = num_samples;
    const int n_dip = n / 2;
    diploid_matrix.resize(static_cast<size_t>(num_sites) * n_dip);
    const int* gmat = genotype_matrix.data();
    int* dmat = diploid_matrix.data();
    for (int s = 0; s < num_sites; s++) {
        const int* g = gmat + static_cast<size_t>(s) * n;
        int* d = dmat + static_cast<size_t>(s) * n_dip;
        for (int i = 0; i < n_dip; i++) {
            d[i] = g[2 * i] + g[2 * i + 1];
        }
    }
    diploid_materialized = true;
}

void ThreadingInstructions::materialize_normalized_haploid() {
    if (standardized_hap_ready) return;
    materialize_genotypes();
    const int n = num_samples;
    standardized_hap.resize(static_cast<size_t>(num_sites) * n);
    const int* gmat = genotype_matrix.data();
    double* zmat = standardized_hap.data();
    for (int s = 0; s < num_sites; s++) {
        const int* g = gmat + static_cast<size_t>(s) * n;
        double* z = zmat + static_cast<size_t>(s) * n;
        int ac = 0;
        for (int i = 0; i < n; i++) ac += g[i];
        const double mu = static_cast<double>(ac) / n;
        const double inv_std = 1.0 / std::sqrt(mu * (1.0 - mu));
        for (int i = 0; i < n; i++) {
            z[i] = (g[i] - mu) * inv_std;
        }
    }
    standardized_hap_ready = true;
}

void ThreadingInstructions::materialize_normalized_diploid() {
    if (standardized_dip_ready) return;
    materialize_diploid();
    const int n = num_samples;
    const int n_dip = n / 2;
    standardized_dip.resize(static_cast<size_t>(num_sites) * n_dip);
    const int* dmat = diploid_matrix.data();
    double* zmat = standardized_dip.data();
    for (int s = 0; s < num_sites; s++) {
        const int* d = dmat + static_cast<size_t>(s) * n_dip;
        double* z = zmat + static_cast<size_t>(s) * n_dip;
        int ac = 0;
        for (int i = 0; i < n_dip; i++) ac += d[i];
        const double mu = static_cast<double>(ac) / n_dip;
        double sample_var = 0.0;
        for (int i = 0; i < n_dip; i++) {
            double diff = d[i] - mu;
            sample_var += diff * diff;
        }
        sample_var /= n_dip;
        const double inv_std = 1.0 / std::sqrt(sample_var);
        for (int i = 0; i < n_dip; i++) {
            z[i] = (d[i] - mu) * inv_std;
        }
    }
    standardized_dip_ready = true;
}

std::vector<double> ThreadingInstructions::left_multiply(const std::vector<double>& x, bool diploid, bool normalize) {
    const int n = num_samples;
    const int n_dip = n / 2;

    if (diploid) {
        if (static_cast<int>(x.size()) != n_dip) {
            std::ostringstream oss;
            oss << "Input vector must have length " << n_dip << ".";
            throw std::runtime_error(oss.str());
        }
    } else {
        if (static_cast<int>(x.size()) != n) {
            std::ostringstream oss;
            oss << "Input vector must have length " << n << ".";
            throw std::runtime_error(oss.str());
        }
    }

    const double* xp = x.data();
    std::vector<double> out(num_sites);
    double* outp = out.data();

    if (normalize && diploid) {
        materialize_normalized_diploid();
        const double* zmat = standardized_dip.data();
        for (int s = 0; s < num_sites; s++) {
            const double* z = zmat + static_cast<size_t>(s) * n_dip;
            double entry = 0.0;
            for (int i = 0; i < n_dip; i++) {
                entry += xp[i] * z[i];
            }
            outp[s] = entry;
        }
    } else if (normalize) {
        materialize_normalized_haploid();
        const double* zmat = standardized_hap.data();
        for (int s = 0; s < num_sites; s++) {
            const double* z = zmat + static_cast<size_t>(s) * n;
            double entry = 0.0;
            for (int i = 0; i < n; i++) {
                entry += xp[i] * z[i];
            }
            outp[s] = entry;
        }
    } else if (diploid) {
        materialize_diploid();
        const int* dmat = diploid_matrix.data();
        for (int s = 0; s < num_sites; s++) {
            const int* d = dmat + static_cast<size_t>(s) * n_dip;
            double entry = 0.0;
            for (int i = 0; i < n_dip; i++) {
                entry += xp[i] * d[i];
            }
            outp[s] = entry;
        }
    } else {
        materialize_genotypes();
        const int* gmat = genotype_matrix.data();
        for (int s = 0; s < num_sites; s++) {
            const int* g = gmat + static_cast<size_t>(s) * n;
            double entry = 0.0;
            for (int i = 0; i < n; i++) {
                entry += xp[i] * g[i];
            }
            outp[s] = entry;
        }
    }
    return out;
}

std::vector<double> ThreadingInstructions::right_multiply(const std::vector<double>& x, bool diploid, bool normalize) {
    if (static_cast<int>(x.size()) != num_sites) {
        std::ostringstream oss;
        oss << "Input vector must have length " << num_sites << ".";
        throw std::runtime_error(oss.str());
    }

    const int n = num_samples;
    const int n_dip = n / 2;
    const int out_size = diploid ? n_dip : n;
    const double* xp = x.data();

    std::vector<double> out(out_size, 0.0);
    double* outp = out.data();

    if (normalize && diploid) {
        materialize_normalized_diploid();
        const double* zmat = standardized_dip.data();
        for (int s = 0; s < num_sites; s++) {
            const double* z = zmat + static_cast<size_t>(s) * n_dip;
            const double w = xp[s];
            for (int i = 0; i < n_dip; i++) {
                outp[i] += w * z[i];
            }
        }
    } else if (normalize) {
        materialize_normalized_haploid();
        const double* zmat = standardized_hap.data();
        for (int s = 0; s < num_sites; s++) {
            const double* z = zmat + static_cast<size_t>(s) * n;
            const double w = xp[s];
            for (int i = 0; i < n; i++) {
                outp[i] += w * z[i];
            }
        }
    } else if (diploid) {
        materialize_diploid();
        const int* dmat = diploid_matrix.data();
        for (int s = 0; s < num_sites; s++) {
            const int* d = dmat + static_cast<size_t>(s) * n_dip;
            const double w = xp[s];
            for (int i = 0; i < n_dip; i++) {
                outp[i] += w * d[i];
            }
        }
    } else {
        materialize_genotypes();
        const int* gmat = genotype_matrix.data();
        for (int s = 0; s < num_sites; s++) {
            const int* g = gmat + static_cast<size_t>(s) * n;
            const double w = xp[s];
            for (int i = 0; i < n; i++) {
                outp[i] += w * g[i];
            }
        }
    }
    return out;
}

void ThreadingInstructions::prepare_tree_multiply() {
    if (tree_ready) return;

    const int n = num_samples;
    const int m = num_sites;

    // ── Phase 1: Build global intervals (unchanged) ─────────────────────

    tree_ref_genome.assign(m, 0);
    for (int idx : instructions[0].mismatches) {
        if (idx >= 0 && idx < m) tree_ref_genome[idx] = 1;
    }

    std::set<int> break_set;
    break_set.insert(0);
    for (int i = 0; i < n; i++) {
        const auto& starts_i = instructions[i].starts;
        const int n_segs = static_cast<int>(starts_i.size());
        for (int k = 1; k < n_segs; k++) {
            auto it = std::lower_bound(positions.begin() + 1, positions.end(), starts_i[k]);
            if (it != positions.end()) {
                break_set.insert(static_cast<int>(it - positions.begin()));
            }
        }
    }

    std::vector<int> breaks(break_set.begin(), break_set.end());
    const int ni = static_cast<int>(breaks.size());
    tree_n_intervals = ni;
    tree_ivl_start.resize(ni);
    tree_ivl_end.resize(ni);
    for (int j = 0; j < ni; j++) {
        tree_ivl_start[j] = breaks[j];
        tree_ivl_end[j] = (j + 1 < ni) ? breaks[j + 1] : m;
    }

    // ── Phase 2: Per-sample segment-to-interval mapping ─────────────────
    // For each sample, map its segments to global intervals via binary search.
    // O(n * S * log(ni)) instead of the previous O(n * ni).

    tree_seg_offset.resize(n);
    tree_seg_first_ivl.clear();
    tree_seg_id.clear();
    for (int i = 0; i < n; i++) {
        tree_seg_offset[i] = static_cast<int>(tree_seg_first_ivl.size());
        const auto& starts_i = instructions[i].starts;
        const int n_segs = static_cast<int>(starts_i.size());
        for (int seg = 0; seg < n_segs; seg++) {
            int first_ivl;
            if (seg == 0) {
                first_ivl = 0;
            } else {
                // Find the site index for this segment start
                auto it = std::lower_bound(positions.begin() + 1, positions.end(),
                                           starts_i[seg]);
                if (it == positions.end()) break;  // beyond all sites
                const int site = static_cast<int>(it - positions.begin());
                // Binary search in breaks to find the interval
                auto jt = std::lower_bound(breaks.begin(), breaks.end(), site);
                first_ivl = static_cast<int>(jt - breaks.begin());
            }
            tree_seg_first_ivl.push_back(first_ivl);
            tree_seg_id.push_back(seg);
        }
    }
    // Sentinel for end-of-last-sample
    tree_seg_first_ivl.push_back(0);
    tree_seg_id.push_back(0);

    // ── Phase 3: Site-to-interval lookup + mismatch signs ──────────────
    // Build O(m) site_to_ivl lookup table, then compute mismatch signs
    // via per-mismatch chain tracing: O(M_total * d * log S) where d = tree
    // depth. Sublinear in n*m when mismatch_rate * d * log S < 1.

    // Build site_to_ivl: O(m)
    site_to_ivl.resize(m);
    for (int j = 0, s = 0; j < ni; j++) {
        while (s < tree_ivl_end[j]) {
            site_to_ivl[s] = j;
            s++;
        }
    }

    // Allocate mismatch sign storage
    tree_mm_offset.resize(n);
    size_t total_mm = 0;
    for (int i = 0; i < n; i++) {
        tree_mm_offset[i] = static_cast<int>(total_mm);
        total_mm += instructions[i].mismatches.size();
    }
    tree_mm_sign.assign(total_mm, 0);

    // Samples 1..n-1: compute mismatch signs via chain tracing
    for (int i = 1; i < n; i++) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const auto& starts_i = instructions[i].starts;

        for (int k = 0; k < n_mm; k++) {
            const int s = mm_data[k];
            const int pos = positions[s];
            int seg = static_cast<int>(
                std::upper_bound(starts_i.begin(), starts_i.end(), pos) - starts_i.begin()) - 1;
            const int target = instructions[i].targets[seg];
            if (target != i) {
                const int g_target = genotype_at_site(target, s);
                tree_mm_sign[tree_mm_offset[i] + k] = static_cast<int8_t>(1 - 2 * g_target);
            }
        }
    }

    // ── Phase 4: Inverted mismatch index, interval changes, self-ref list ──

    // 4a: Inverted mismatch index — group non-self-ref mismatches by interval.
    // Count per interval first.
    std::vector<int> inv_count(ni, 0);
    for (int i = 1; i < n; i++) {
        const int off = tree_mm_offset[i];
        const int n_mm_i = static_cast<int>(instructions[i].mismatches.size());
        for (int k = 0; k < n_mm_i; k++) {
            if (tree_mm_sign[off + k] != 0) {
                inv_count[site_to_ivl[instructions[i].mismatches[k]]]++;
            }
        }
    }
    inv_mm_offset.resize(ni + 1);
    inv_mm_offset[0] = 0;
    for (int j = 0; j < ni; j++) inv_mm_offset[j + 1] = inv_mm_offset[j] + inv_count[j];
    const int total_inv = inv_mm_offset[ni];
    inv_mm_sample.resize(total_inv);
    inv_mm_site.resize(total_inv);
    inv_mm_sign.resize(total_inv);
    std::fill(inv_count.begin(), inv_count.end(), 0);
    for (int i = 1; i < n; i++) {
        const int off = tree_mm_offset[i];
        const int n_mm_i = static_cast<int>(instructions[i].mismatches.size());
        const int* mm_data = instructions[i].mismatches.data();
        for (int k = 0; k < n_mm_i; k++) {
            if (tree_mm_sign[off + k] != 0) {
                const int j = site_to_ivl[mm_data[k]];
                const int idx = inv_mm_offset[j] + inv_count[j]++;
                inv_mm_sample[idx] = i;
                inv_mm_site[idx] = mm_data[k];
                inv_mm_sign[idx] = tree_mm_sign[off + k];
            }
        }
    }

    // 4b: Interval change list — which samples change target at each boundary.
    std::vector<int> chg_count(ni, 0);
    for (int i = 0; i < n; i++) {
        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped = next_off - seg_off;
        for (int sk = 1; sk < n_mapped; sk++) {
            const int boundary_ivl = tree_seg_first_ivl[seg_off + sk];
            if (boundary_ivl > 0) chg_count[boundary_ivl]++;
        }
    }
    ivl_change_offset.resize(ni + 1);
    ivl_change_offset[0] = 0;
    for (int j = 0; j < ni; j++) ivl_change_offset[j + 1] = ivl_change_offset[j] + chg_count[j];
    const int total_changes = ivl_change_offset[ni];
    ivl_change_sample.resize(total_changes);
    ivl_change_old_tgt.resize(total_changes);
    ivl_change_new_tgt.resize(total_changes);
    std::fill(chg_count.begin(), chg_count.end(), 0);
    for (int i = 0; i < n; i++) {
        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped = next_off - seg_off;
        for (int sk = 1; sk < n_mapped; sk++) {
            const int boundary_ivl = tree_seg_first_ivl[seg_off + sk];
            if (boundary_ivl > 0) {
                const int old_seg = tree_seg_id[seg_off + sk - 1];
                const int new_seg = tree_seg_id[seg_off + sk];
                const int idx = ivl_change_offset[boundary_ivl] + chg_count[boundary_ivl]++;
                ivl_change_sample[idx] = i;
                ivl_change_old_tgt[idx] = instructions[i].targets[old_seg];
                ivl_change_new_tgt[idx] = instructions[i].targets[new_seg];
            }
        }
    }
    // Reverse each boundary group to get descending sample order (for bottom-up detach).
    // Data was filled in increasing order since the outer loop iterates i = 0..n-1.
    for (int j = 0; j < ni; j++) {
        int a = ivl_change_offset[j], b = ivl_change_offset[j + 1] - 1;
        while (a < b) {
            std::swap(ivl_change_sample[a], ivl_change_sample[b]);
            std::swap(ivl_change_old_tgt[a], ivl_change_old_tgt[b]);
            std::swap(ivl_change_new_tgt[a], ivl_change_new_tgt[b]);
            a++; b--;
        }
    }

    // 4c: Self-ref samples per interval (sample 0 excluded — handled via ref term)
    std::vector<int> sr_count(ni, 0);
    for (int i = 1; i < n; i++) {
        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped = next_off - seg_off;
        for (int sk = 0; sk < n_mapped; sk++) {
            const int seg = tree_seg_id[seg_off + sk];
            if (instructions[i].targets[seg] == i) {
                const int fi = tree_seg_first_ivl[seg_off + sk];
                const int li = (sk + 1 < n_mapped) ? tree_seg_first_ivl[seg_off + sk + 1] : ni;
                for (int j = fi; j < li; j++) sr_count[j]++;
            }
        }
    }
    self_ref_offset.resize(ni + 1);
    self_ref_offset[0] = 0;
    for (int j = 0; j < ni; j++) self_ref_offset[j + 1] = self_ref_offset[j] + sr_count[j];
    const int total_sr = self_ref_offset[ni];
    self_ref_sample.resize(total_sr);
    std::fill(sr_count.begin(), sr_count.end(), 0);
    for (int i = 1; i < n; i++) {
        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped = next_off - seg_off;
        for (int sk = 0; sk < n_mapped; sk++) {
            const int seg = tree_seg_id[seg_off + sk];
            if (instructions[i].targets[seg] == i) {
                const int fi = tree_seg_first_ivl[seg_off + sk];
                const int li = (sk + 1 < n_mapped) ? tree_seg_first_ivl[seg_off + sk + 1] : ni;
                for (int j = fi; j < li; j++) {
                    self_ref_sample[self_ref_offset[j] + sr_count[j]++] = i;
                }
            }
        }
    }

    // 4d: Precompute ancestor chains for left multiply W updates.
    // For each change at boundary b:
    //   detach chain: trace from old_tgt up to root using tree at interval b-1
    //   attach chain: trace from new_tgt up to root using tree at interval b
    detach_chain_offset.resize(total_changes + 1);
    attach_chain_offset.resize(total_changes + 1);
    detach_chain_data.clear();
    attach_chain_data.clear();
    detach_chain_offset[0] = 0;
    attach_chain_offset[0] = 0;

    for (int b = 1; b < ni; b++) {
        for (int c = ivl_change_offset[b]; c < ivl_change_offset[b + 1]; c++) {
            const int sample_c = ivl_change_sample[c];
            const int old_tgt = ivl_change_old_tgt[c];
            const int new_tgt = ivl_change_new_tgt[c];

            // Detach chain: old tree (interval b-1)
            if (old_tgt != sample_c) {
                int cur = old_tgt;
                while (true) {
                    detach_chain_data.push_back(cur);
                    if (cur == 0) break;
                    const int next = target_at_interval(cur, b - 1);
                    if (next == cur) break;
                    cur = next;
                }
            }
            detach_chain_offset[c + 1] = static_cast<int>(detach_chain_data.size());

            // Attach chain: new tree (interval b)
            if (new_tgt != sample_c) {
                int cur = new_tgt;
                while (true) {
                    attach_chain_data.push_back(cur);
                    if (cur == 0) break;
                    const int next = target_at_interval(cur, b);
                    if (next == cur) break;
                    cur = next;
                }
            }
            attach_chain_offset[c + 1] = static_cast<int>(attach_chain_data.size());
        }
    }

    tree_ready = true;
}

int ThreadingInstructions::target_at_interval(int sample, int interval) const {
    const int seg_off = tree_seg_offset[sample];
    const int next_off = (sample + 1 < num_samples) ? tree_seg_offset[sample + 1]
                                                     : static_cast<int>(tree_seg_first_ivl.size()) - 1;
    // Binary search: find last mapped segment whose first_ivl <= interval
    const int* base = tree_seg_first_ivl.data() + seg_off;
    const int n_mapped = next_off - seg_off;
    int lo = 0, hi = n_mapped - 1;
    while (lo < hi) {
        int mid = lo + (hi - lo + 1) / 2;
        if (base[mid] <= interval) lo = mid;
        else hi = mid - 1;
    }
    return instructions[sample].targets[tree_seg_id[seg_off + lo]];
}

int ThreadingInstructions::genotype_at_site(int sample, int site) const {
    int parity = 0;
    int cur = sample;
    while (cur > 0) {
        const auto& instr = instructions[cur];
        int seg = static_cast<int>(
            std::upper_bound(instr.starts.begin(), instr.starts.end(),
                             positions[site]) - instr.starts.begin()) - 1;
        if (std::binary_search(instr.mismatches.begin(),
                               instr.mismatches.end(), site)) {
            parity ^= 1;
        }
        int next = instr.targets[seg];
        if (next == cur) break;  // self-ref: stop tracing
        cur = next;
    }
    return tree_ref_genome[site] ^ parity;
}

std::vector<double> ThreadingInstructions::right_multiply_tree(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_sites) {
        std::ostringstream oss;
        oss << "Input vector must have length " << num_sites << ".";
        throw std::runtime_error(oss.str());
    }

    prepare_tree_multiply();

    const int n = num_samples;
    const int ni = tree_n_intervals;
    const double* xp = x.data();
    const int* ref = tree_ref_genome.data();
    const int* ivl_s = tree_ivl_start.data();
    const int* ivl_e = tree_ivl_end.data();

    std::vector<double> out(n, 0.0);

    // Region-chunked parallel sweep: each thread handles a contiguous
    // range of intervals with its own tree state and local output.
    #pragma omp parallel
    {
        const int T = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        const int chunk_size = (ni + T - 1) / T;
        const int j_start = std::min(tid * chunk_size, ni);
        const int j_end = std::min(j_start + chunk_size, ni);

        if (j_start < j_end) {
            // Thread-local state
            std::vector<int> cur_target(n);
            std::vector<double> delta(n, 0.0);
            std::vector<double> mm_corr(n, 0.0);
            std::vector<double> local_out(n, 0.0);
            double total_ref = 0.0;

            // Initialize cur_target for this chunk's starting interval
            cur_target[0] = 0;
            for (int i = 1; i < n; i++) {
                cur_target[i] = target_at_interval(i, j_start);
            }

            // Sweep intervals [j_start, j_end)
            for (int j = j_start; j < j_end; j++) {
                // Apply target changes (skip for j_start: already in cur_target)
                if (j > j_start) {
                    for (int c = ivl_change_offset[j]; c < ivl_change_offset[j + 1]; c++) {
                        cur_target[ivl_change_sample[c]] = ivl_change_new_tgt[c];
                    }
                }

                // Ref contribution
                double ref_sum = 0.0;
                for (int s = ivl_s[j]; s < ivl_e[j]; s++) {
                    ref_sum += xp[s] * ref[s];
                }
                total_ref += ref_sum;

                // Skip if no corrections in this interval
                const int mm_a = inv_mm_offset[j];
                const int mm_b = inv_mm_offset[j + 1];
                if (mm_a == mm_b) continue;

                // Scatter mismatch corrections
                for (int c = mm_a; c < mm_b; c++) {
                    mm_corr[inv_mm_sample[c]] += xp[inv_mm_site[c]] * inv_mm_sign[c];
                }

                // Propagate corrections from root to leaves.
                // delta[i] is fully overwritten (target[i]<i ensures deps ready),
                // so no fill between intervals is needed.
                delta[0] = 0.0;
                for (int i = 1; i < n; i++) {
                    const double d = delta[cur_target[i]] + mm_corr[i];
                    delta[i] = d;
                    local_out[i] += d;
                }

                // Clear mm_corr sparsely
                for (int c = mm_a; c < mm_b; c++) {
                    mm_corr[inv_mm_sample[c]] = 0.0;
                }
            }

            // Add deferred ref contribution
            for (int i = 0; i < n; i++) local_out[i] += total_ref;

            // Reduce into shared output
            #pragma omp critical
            {
                for (int i = 0; i < n; i++) out[i] += local_out[i];
            }
        }
    }

    return out;
}

std::vector<double> ThreadingInstructions::left_multiply_tree(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_samples) {
        std::ostringstream oss;
        oss << "Input vector must have length " << num_samples << ".";
        throw std::runtime_error(oss.str());
    }

    prepare_tree_multiply();

    const int n = num_samples;
    const int m = num_sites;
    const int ni = tree_n_intervals;
    const double* xp = x.data();

    // Sample 0's ref=1 sites (sorted) for fast ref-term emission
    const int* ref0_mm = instructions[0].mismatches.data();
    const int n_ref0 = static_cast<int>(instructions[0].mismatches.size());

    std::vector<double> out(m, 0.0);

    // Region-chunked parallel sweep: each thread builds its own W vector
    // for its starting interval and sweeps its chunk. Threads write to
    // non-overlapping site ranges in out[], so no reduction needed.
    #pragma omp parallel
    {
        const int T = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        const int chunk_size = (ni + T - 1) / T;
        const int j_start = std::min(tid * chunk_size, ni);
        const int j_end = std::min(j_start + chunk_size, ni);

        if (j_start < j_end) {
            // Build W for this chunk's starting interval
            std::vector<double> W(n);
            for (int i = 0; i < n; i++) W[i] = xp[i];
            for (int i = n - 1; i >= 1; i--) {
                const int t = target_at_interval(i, j_start);
                if (t != i) W[t] += W[i];
            }

            // Preallocated buffer for W-update snapshots
            std::vector<double> saved_delta;
            double* outp = out.data();

            // Sweep intervals [j_start, j_end)
            for (int j = j_start; j < j_end; j++) {
                const int site_a = tree_ivl_start[j];
                const int site_b = tree_ivl_end[j];

                // A. Ref contribution: out[s] += W[0] for ref[s]==1 sites
                {
                    const int* lo = std::lower_bound(ref0_mm, ref0_mm + n_ref0, site_a);
                    const int* hi = std::lower_bound(lo, ref0_mm + n_ref0, site_b);
                    const double w0 = W[0];
                    for (const int* it = lo; it != hi; ++it) {
                        outp[*it] += w0;
                    }
                }

                // B. Non-self-ref mismatch contributions
                {
                    const int mm_a = inv_mm_offset[j];
                    const int mm_b = inv_mm_offset[j + 1];
                    for (int k = mm_a; k < mm_b; k++) {
                        outp[inv_mm_site[k]] += static_cast<double>(inv_mm_sign[k]) * W[inv_mm_sample[k]];
                    }
                }

                // C. Self-ref carry-run (samples > 0 only)
                {
                    const int sr_a = self_ref_offset[j];
                    const int sr_b = self_ref_offset[j + 1];
                    for (int k = sr_a; k < sr_b; k++) {
                        const int i = self_ref_sample[k];
                        const double wi = W[i];
                        const int* mm_data = instructions[i].mismatches.data();
                        const int n_mm_i = static_cast<int>(instructions[i].mismatches.size());

                        int carry = (j == 0) ? 0
                                    : genotype_at_site(i, tree_ivl_end[j - 1] - 1);
                        const int* lo = std::lower_bound(mm_data, mm_data + n_mm_i, site_a);
                        const int* hi_mm = std::lower_bound(mm_data, mm_data + n_mm_i, site_b);
                        int prev = site_a;

                        for (const int* it = lo; it != hi_mm; ++it) {
                            const int mm_site = *it;
                            if (carry && mm_site > prev) {
                                double* __restrict__ op = outp + prev;
                                const int len = mm_site - prev;
                                for (int s = 0; s < len; s++) op[s] += wi;
                            }
                            carry ^= 1;
                            outp[mm_site] += wi * carry;
                            prev = mm_site + 1;
                        }
                        if (carry && site_b > prev) {
                            double* __restrict__ op = outp + prev;
                            const int len = site_b - prev;
                            for (int s = 0; s < len; s++) op[s] += wi;
                        }
                    }
                }

                // D. Update W at boundary j → j+1 using precomputed ancestor chains.
                if (j + 1 < j_end) {
                    const int chg_a = ivl_change_offset[j + 1];
                    const int chg_b = ivl_change_offset[j + 2];
                    if (chg_a < chg_b) {
                        const int n_chg = chg_b - chg_a;

                        // Pass 1: Detach via precomputed chains (descending sample order)
                        for (int k = chg_a; k < chg_b; k++) {
                            const int dc_a = detach_chain_offset[k];
                            const int dc_b = detach_chain_offset[k + 1];
                            if (dc_a == dc_b) continue;
                            const double d = W[ivl_change_sample[k]];
                            for (int p = dc_a; p < dc_b; p++) {
                                W[detach_chain_data[p]] -= d;
                            }
                        }

                        // Snapshot deltas after detach, before attach
                        saved_delta.resize(n_chg);
                        for (int k = 0; k < n_chg; k++) {
                            saved_delta[k] = W[ivl_change_sample[chg_a + k]];
                        }

                        // Pass 2: Attach via precomputed chains
                        for (int k = chg_a; k < chg_b; k++) {
                            const int ac_a = attach_chain_offset[k];
                            const int ac_b = attach_chain_offset[k + 1];
                            if (ac_a == ac_b) continue;
                            const double d = saved_delta[k - chg_a];
                            for (int p = ac_a; p < ac_b; p++) {
                                W[attach_chain_data[p]] += d;
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}


// ═══════════════════════════════════════════════════════════════════════════
// Batch tree multiply: process k vectors in a single tree traversal.
// Layout: row-major flat arrays, X[row * k + col].
// ═══════════════════════════════════════════════════════════════════════════

std::vector<double> ThreadingInstructions::right_multiply_tree_batch(
        const std::vector<double>& x_flat, int k) {
    if (static_cast<int>(x_flat.size()) != num_sites * k)
        throw std::runtime_error("Input must have length num_sites * k");
    if (k == 1) return right_multiply_tree(x_flat);

    // Sequential per-column calls: each call uses interval-parallel OMP
    // internally with O(n) working set per thread.  This keeps the delta
    // array small enough to fit in L2 cache, which is critical for the
    // random-access pattern delta[cur_target[i]].  A k-wide interleaved
    // approach was tested but regresses at large n because the expanded
    // n*k delta array exceeds L2, increasing cache miss cost.
    const int n = num_samples;
    const int m = num_sites;
    const double* xp = x_flat.data();

    std::vector<double> out(static_cast<size_t>(n) * k);

    for (int col = 0; col < k; col++) {
        std::vector<double> x_col(m);
        for (int s = 0; s < m; s++)
            x_col[s] = xp[static_cast<size_t>(s) * k + col];

        auto result = right_multiply_tree(x_col);

        for (int i = 0; i < n; i++)
            out[static_cast<size_t>(i) * k + col] = result[i];
    }

    return out;
}


std::vector<double> ThreadingInstructions::left_multiply_tree_batch(
        const std::vector<double>& x_flat, int k) {
    if (static_cast<int>(x_flat.size()) != num_samples * k)
        throw std::runtime_error("Input must have length num_samples * k");
    if (k == 1) return left_multiply_tree(x_flat);

    prepare_tree_multiply();
    const int n = num_samples, m = num_sites, ni = tree_n_intervals;
    const double* xp = x_flat.data();

    const int* ref0_mm = instructions[0].mismatches.data();
    const int n_ref0 = static_cast<int>(instructions[0].mismatches.size());

    std::vector<double> out(static_cast<size_t>(m) * k, 0.0);

    // Region-chunked parallel sweep with k-wide W vectors.  Each thread
    // builds its own W for its starting interval and sweeps its chunk.
    // Threads write to non-overlapping site ranges in out[].
    #pragma omp parallel
    {
        const int T = omp_get_num_threads();
        const int tid = omp_get_thread_num();
        const int chunk_size = (ni + T - 1) / T;
        const int j_start = std::min(tid * chunk_size, ni);
        const int j_end = std::min(j_start + chunk_size, ni);

        if (j_start < j_end) {
            // Build k-wide W for this chunk's starting interval
            const size_t nk = static_cast<size_t>(n) * k;
            std::vector<double> W(nk);
            for (size_t idx = 0; idx < nk; idx++) W[idx] = xp[idx];
            for (int i = n - 1; i >= 1; i--) {
                const int t = target_at_interval(i, j_start);
                if (t != i) {
                    const double* wi = &W[static_cast<size_t>(i) * k];
                    double* wt = &W[static_cast<size_t>(t) * k];
                    for (int c = 0; c < k; c++) wt[c] += wi[c];
                }
            }

            std::vector<double> saved_delta;
            double* outp = out.data();

            for (int j = j_start; j < j_end; j++) {
                const int site_a = tree_ivl_start[j];
                const int site_b = tree_ivl_end[j];

                // A. Ref contribution
                {
                    const int* lo = std::lower_bound(ref0_mm, ref0_mm + n_ref0, site_a);
                    const int* hi = std::lower_bound(lo, ref0_mm + n_ref0, site_b);
                    const double* w0 = &W[0];
                    for (const int* it = lo; it != hi; ++it) {
                        double* od = &outp[static_cast<size_t>(*it) * k];
                        for (int c = 0; c < k; c++) od[c] += w0[c];
                    }
                }

                // B. Non-self-ref mismatch contributions
                {
                    const int mm_a = inv_mm_offset[j];
                    const int mm_b = inv_mm_offset[j + 1];
                    for (int kk = mm_a; kk < mm_b; kk++) {
                        const double sv = static_cast<double>(inv_mm_sign[kk]);
                        const double* wi = &W[static_cast<size_t>(inv_mm_sample[kk]) * k];
                        double* od = &outp[static_cast<size_t>(inv_mm_site[kk]) * k];
                        for (int c = 0; c < k; c++) od[c] += wi[c] * sv;
                    }
                }

                // C. Self-ref carry-run
                {
                    const int sr_a = self_ref_offset[j];
                    const int sr_b = self_ref_offset[j + 1];
                    for (int kk = sr_a; kk < sr_b; kk++) {
                        const int i = self_ref_sample[kk];
                        const double* wi = &W[static_cast<size_t>(i) * k];
                        const int* mm_data = instructions[i].mismatches.data();
                        const int n_mm_i = static_cast<int>(instructions[i].mismatches.size());

                        int carry = (j == 0) ? 0
                                    : genotype_at_site(i, tree_ivl_end[j - 1] - 1);
                        const int* lo = std::lower_bound(mm_data, mm_data + n_mm_i, site_a);
                        const int* hi_mm = std::lower_bound(mm_data, mm_data + n_mm_i, site_b);
                        int prev = site_a;

                        for (const int* it = lo; it != hi_mm; ++it) {
                            const int ms = *it;
                            if (carry && ms > prev) {
                                for (int s = prev; s < ms; s++) {
                                    double* od = &outp[static_cast<size_t>(s) * k];
                                    for (int c = 0; c < k; c++) od[c] += wi[c];
                                }
                            }
                            carry ^= 1;
                            if (carry) {
                                double* od = &outp[static_cast<size_t>(ms) * k];
                                for (int c = 0; c < k; c++) od[c] += wi[c];
                            }
                            prev = ms + 1;
                        }
                        if (carry && site_b > prev) {
                            for (int s = prev; s < site_b; s++) {
                                double* od = &outp[static_cast<size_t>(s) * k];
                                for (int c = 0; c < k; c++) od[c] += wi[c];
                            }
                        }
                    }
                }

                // D. Update W at boundary j → j+1
                if (j + 1 < j_end) {
                    const int chg_a = ivl_change_offset[j + 1];
                    const int chg_b = ivl_change_offset[j + 2];
                    if (chg_a < chg_b) {
                        const int n_chg = chg_b - chg_a;

                        // Pass 1: Detach via precomputed chains
                        for (int kk = chg_a; kk < chg_b; kk++) {
                            const int dc_a = detach_chain_offset[kk];
                            const int dc_b = detach_chain_offset[kk + 1];
                            if (dc_a == dc_b) continue;
                            const double* d = &W[static_cast<size_t>(ivl_change_sample[kk]) * k];
                            for (int p = dc_a; p < dc_b; p++) {
                                double* wc = &W[static_cast<size_t>(detach_chain_data[p]) * k];
                                for (int c = 0; c < k; c++) wc[c] -= d[c];
                            }
                        }

                        // Snapshot deltas after detach
                        saved_delta.resize(static_cast<size_t>(n_chg) * k);
                        for (int kk = 0; kk < n_chg; kk++) {
                            const double* wi = &W[static_cast<size_t>(ivl_change_sample[chg_a + kk]) * k];
                            double* sd = &saved_delta[static_cast<size_t>(kk) * k];
                            for (int c = 0; c < k; c++) sd[c] = wi[c];
                        }

                        // Pass 2: Attach via precomputed chains
                        for (int kk = chg_a; kk < chg_b; kk++) {
                            const int ac_a = attach_chain_offset[kk];
                            const int ac_b = attach_chain_offset[kk + 1];
                            if (ac_a == ac_b) continue;
                            const double* d = &saved_delta[static_cast<size_t>(kk - chg_a) * k];
                            for (int p = ac_a; p < ac_b; p++) {
                                double* wc = &W[static_cast<size_t>(attach_chain_data[p]) * k];
                                for (int c = 0; c < k; c++) wc[c] += d[c];
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}



// ---------------------------------------------------------------------------
// ARG traversal: visit_clades / visit_branches
// ---------------------------------------------------------------------------
//
// Reconstructs the implicit coalescent tree from threading instructions at each
// interval (contiguous region where no segment boundary changes the tree) and
// enumerates clades/branches via union-find.

namespace {

struct SimpleUF {
    std::vector<int> par;
    std::vector<int> rnk;
    std::vector<std::vector<int>> members;
    std::vector<double> height;  // height of each root's subtree

    explicit SimpleUF(int n) : par(n), rnk(n, 0), members(n), height(n, 0.0) {
        for (int i = 0; i < n; i++) {
            par[i] = i;
            members[i] = {i};
        }
    }

    int find(int x) {
        while (par[x] != par[par[x]]) par[x] = par[par[x]];
        while (par[x] != x) x = par[x];
        return x;
    }

    // Merge sets containing a and b.  Returns root of merged set.
    int unite(int a, int b) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return ra;
        if (rnk[ra] < rnk[rb]) std::swap(ra, rb);
        par[rb] = ra;
        if (rnk[ra] == rnk[rb]) rnk[ra]++;
        members[ra].insert(members[ra].end(), members[rb].begin(), members[rb].end());
        members[rb].clear();
        return ra;
    }

    void reset(int n) {
        for (int i = 0; i < n; i++) {
            par[i] = i;
            rnk[i] = 0;
            members[i].clear();
            members[i].push_back(i);
            height[i] = 0.0;
        }
    }
};

struct Edge {
    int child;
    int parent;
    double height;
};

// Find segment index for a sample at a given bp position via binary search.
// instructions[i].starts contains bp positions.
inline int segment_at_bp(const ThreadingInstruction& instr, int bp) {
    if (instr.num_segments == 0) return -1;
    auto it = std::upper_bound(instr.starts.begin(), instr.starts.end(), bp);
    return static_cast<int>(it - instr.starts.begin()) - 1;
}

} // anonymous namespace


void ThreadingInstructions::visit_clades(
    std::function<void(const std::vector<int>&, double, int, int)> callback) const
{
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return;

    // Collect all segment boundary bp positions across all samples
    std::set<int> boundary_set;
    boundary_set.insert(start);
    for (int i = 0; i < n; i++) {
        for (int s : instructions[i].starts) {
            if (s > start && s < end)
                boundary_set.insert(s);
        }
    }
    std::vector<int> boundaries(boundary_set.begin(), boundary_set.end());

    // Emit leaf clades once, spanning the full range
    for (int i = 0; i < n; i++) {
        callback({i}, 0.0, start, end);
    }

    SimpleUF uf(n);
    std::vector<Edge> edges;
    edges.reserve(n);

    // Track previous interval's (target, tmrca) to skip re-processing identical trees
    std::vector<int> prev_target(n, -2);
    std::vector<double> prev_tmrca(n, -1.0);
    std::vector<int> cur_target(n);
    std::vector<double> cur_tmrca(n);

    // Pending clades from the current identical-tree run
    struct PendingClade { std::vector<int> desc; double height; int start_bp; };
    std::vector<PendingClade> pending;

    auto flush_pending = [&](int end_bp) {
        for (auto& pc : pending)
            callback(pc.desc, pc.height, pc.start_bp, end_bp);
        pending.clear();
    };

    for (int b = 0; b < static_cast<int>(boundaries.size()); b++) {
        int bp_pos = boundaries[b];
        int bp_next = (b + 1 < static_cast<int>(boundaries.size()))
                        ? boundaries[b + 1] : end;

        // Determine (target, tmrca) for each sample at this bp position
        cur_target[0] = -1;
        cur_tmrca[0] = 0.0;
        for (int i = 1; i < n; i++) {
            int seg = segment_at_bp(instructions[i], bp_pos);
            if (seg < 0) { cur_target[i] = 0; cur_tmrca[i] = 0.0; continue; }
            cur_target[i] = instructions[i].targets[seg];
            cur_tmrca[i] = instructions[i].tmrcas[seg];
        }

        // Check if tree topology changed from previous interval
        bool same = (b > 0);
        if (same) {
            for (int i = 1; i < n; i++) {
                if (cur_target[i] != prev_target[i] || cur_tmrca[i] != prev_tmrca[i]) {
                    same = false;
                    break;
                }
            }
        }

        if (same) continue;

        // Tree changed — flush previous run's clades
        if (!pending.empty()) flush_pending(bp_pos);

        // Save current state
        std::swap(prev_target, cur_target);
        std::swap(prev_tmrca, cur_tmrca);

        // Build edges from this interval's threading state
        edges.clear();
        for (int i = 1; i < n; i++) {
            if (prev_target[i] >= 0) {
                edges.push_back({i, prev_target[i], prev_tmrca[i]});
            }
        }
        std::sort(edges.begin(), edges.end(),
                  [](const Edge& a, const Edge& b) { return a.height < b.height; });

        // Union-find: build tree bottom-up, emit clades
        uf.reset(n);
        for (const auto& e : edges) {
            int rc = uf.find(e.child);
            int rp = uf.find(e.parent);
            if (rc != rp) {
                int root = uf.unite(rc, rp);
                uf.height[root] = e.height;
                auto& m = uf.members[root];
                std::sort(m.begin(), m.end());
                pending.push_back({m, e.height, bp_pos});
            }
        }
    }

    // Flush final run
    flush_pending(end);
}


void ThreadingInstructions::visit_branches(
    std::function<void(const std::vector<int>&, double, double, int, int)> callback) const
{
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return;

    // Collect all segment boundary bp positions
    std::set<int> boundary_set;
    boundary_set.insert(start);
    for (int i = 0; i < n; i++) {
        for (int s : instructions[i].starts) {
            if (s > start && s < end)
                boundary_set.insert(s);
        }
    }
    std::vector<int> boundaries(boundary_set.begin(), boundary_set.end());

    SimpleUF uf(n);
    std::vector<Edge> edges;
    edges.reserve(n);

    std::vector<int> prev_target(n, -2);
    std::vector<double> prev_tmrca(n, -1.0);
    std::vector<int> cur_target(n);
    std::vector<double> cur_tmrca(n);

    struct PendingBranch { std::vector<int> desc; double child_h; double parent_h; int start_bp; };
    std::vector<PendingBranch> pending;

    auto flush_pending = [&](int end_bp) {
        for (auto& pb : pending)
            callback(pb.desc, pb.child_h, pb.parent_h, pb.start_bp, end_bp);
        pending.clear();
    };

    for (int b = 0; b < static_cast<int>(boundaries.size()); b++) {
        int bp_pos = boundaries[b];

        cur_target[0] = -1;
        cur_tmrca[0] = 0.0;
        for (int i = 1; i < n; i++) {
            int seg = segment_at_bp(instructions[i], bp_pos);
            if (seg < 0) { cur_target[i] = 0; cur_tmrca[i] = 0.0; continue; }
            cur_target[i] = instructions[i].targets[seg];
            cur_tmrca[i] = instructions[i].tmrcas[seg];
        }

        bool same = (b > 0);
        if (same) {
            for (int i = 1; i < n; i++) {
                if (cur_target[i] != prev_target[i] || cur_tmrca[i] != prev_tmrca[i]) {
                    same = false;
                    break;
                }
            }
        }

        if (same) continue;

        if (!pending.empty()) flush_pending(bp_pos);

        std::swap(prev_target, cur_target);
        std::swap(prev_tmrca, cur_tmrca);

        edges.clear();
        for (int i = 1; i < n; i++) {
            if (prev_target[i] >= 0)
                edges.push_back({i, prev_target[i], prev_tmrca[i]});
        }
        std::sort(edges.begin(), edges.end(),
                  [](const Edge& a, const Edge& b) { return a.height < b.height; });

        uf.reset(n);
        for (const auto& e : edges) {
            int rc = uf.find(e.child);
            int rp = uf.find(e.parent);
            if (rc != rp) {
                if (uf.height[rc] < e.height) {
                    auto desc_c = uf.members[rc];
                    std::sort(desc_c.begin(), desc_c.end());
                    pending.push_back({desc_c, uf.height[rc], e.height, bp_pos});
                }

                if (uf.height[rp] < e.height) {
                    auto desc_p = uf.members[rp];
                    std::sort(desc_p.begin(), desc_p.end());
                    pending.push_back({desc_p, uf.height[rp], e.height, bp_pos});
                }

                int root = uf.unite(rc, rp);
                uf.height[root] = e.height;
            }
        }
    }

    flush_pending(end);
}


double ThreadingInstructions::total_volume() const {
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return 0.0;

    // Collect segment boundary bp positions
    std::set<int> boundary_set;
    boundary_set.insert(start);
    for (int i = 0; i < n; i++) {
        for (int s : instructions[i].starts) {
            if (s > start && s < end)
                boundary_set.insert(s);
        }
    }
    std::vector<int> boundaries(boundary_set.begin(), boundary_set.end());

    // Lightweight union-find: only tracks component heights, no member lists
    std::vector<int> parent_uf(n);
    std::vector<int> rank_uf(n, 0);
    std::vector<double> comp_height(n, 0.0);

    auto find = [&](int x) {
        while (parent_uf[x] != x) { parent_uf[x] = parent_uf[parent_uf[x]]; x = parent_uf[x]; }
        return x;
    };
    auto unite = [&](int a, int b) -> int {
        if (rank_uf[a] < rank_uf[b]) std::swap(a, b);
        parent_uf[b] = a;
        if (rank_uf[a] == rank_uf[b]) rank_uf[a]++;
        return a;
    };

    struct LightEdge { int child; int parent; double height; };
    std::vector<LightEdge> edges;
    edges.reserve(n);

    std::vector<int> prev_target(n, -2);
    std::vector<double> prev_tmrca(n, -1.0);

    double vol = 0.0;
    double prev_vol = 0.0;

    for (int b = 0; b < static_cast<int>(boundaries.size()); b++) {
        int bp_pos = boundaries[b];
        int bp_next = (b + 1 < static_cast<int>(boundaries.size()))
                        ? boundaries[b + 1] : end;
        double span = static_cast<double>(bp_next - bp_pos);

        // Check if topology changed
        bool same = (b > 0);
        edges.clear();
        for (int i = 1; i < n; i++) {
            int seg = segment_at_bp(instructions[i], bp_pos);
            int tgt = (seg < 0) ? 0 : instructions[i].targets[seg];
            double tmr = (seg < 0) ? 0.0 : instructions[i].tmrcas[seg];
            if (same && (tgt != prev_target[i] || tmr != prev_tmrca[i]))
                same = false;
            prev_target[i] = tgt;
            prev_tmrca[i] = tmr;
            if (tgt >= 0)
                edges.push_back({i, tgt, tmr});
        }

        if (same) {
            vol += prev_vol * span;
            continue;
        }

        // Build tree, sum branch lengths
        for (int i = 0; i < n; i++) { parent_uf[i] = i; rank_uf[i] = 0; comp_height[i] = 0.0; }
        std::sort(edges.begin(), edges.end(),
                  [](const LightEdge& a, const LightEdge& b) { return a.height < b.height; });

        double interval_vol = 0.0;
        for (const auto& e : edges) {
            int rc = find(e.child);
            int rp = find(e.parent);
            if (rc != rp) {
                interval_vol += (e.height - comp_height[rc]) + (e.height - comp_height[rp]);
                int root = unite(rc, rp);
                comp_height[root] = e.height;
            }
        }
        prev_vol = interval_vol;
        vol += interval_vol * span;
    }
    return vol;
}


std::vector<double> ThreadingInstructions::allele_frequency_spectrum() const {
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return std::vector<double>(n + 1, 0.0);

    std::set<int> boundary_set;
    boundary_set.insert(start);
    for (int i = 0; i < n; i++)
        for (int s : instructions[i].starts)
            if (s > start && s < end)
                boundary_set.insert(s);
    std::vector<int> boundaries(boundary_set.begin(), boundary_set.end());

    // Union-find with component heights and sizes (no member lists)
    std::vector<int> uf_parent(n), uf_rank(n, 0), uf_size(n, 1);
    std::vector<double> uf_height(n, 0.0);

    auto find = [&](int x) {
        while (uf_parent[x] != x) { uf_parent[x] = uf_parent[uf_parent[x]]; x = uf_parent[x]; }
        return x;
    };
    auto unite = [&](int a, int b) -> int {
        if (uf_rank[a] < uf_rank[b]) std::swap(a, b);
        uf_parent[b] = a;
        uf_size[a] += uf_size[b];
        if (uf_rank[a] == uf_rank[b]) uf_rank[a]++;
        return a;
    };

    struct LightEdge { int child; int parent; double height; };
    std::vector<LightEdge> edges;
    edges.reserve(n);

    std::vector<int> prev_target(n, -2);
    std::vector<double> prev_tmrca(n, -1.0);

    std::vector<double> afs(n + 1, 0.0);
    // Cache per-interval AFS contribution (volume per bin, unscaled by span)
    std::vector<double> prev_afs_contrib(n + 1, 0.0);

    for (int b = 0; b < static_cast<int>(boundaries.size()); b++) {
        int bp_pos = boundaries[b];
        int bp_next = (b + 1 < static_cast<int>(boundaries.size()))
                        ? boundaries[b + 1] : end;
        double span = static_cast<double>(bp_next - bp_pos);

        bool same = (b > 0);
        edges.clear();
        for (int i = 1; i < n; i++) {
            int seg = segment_at_bp(instructions[i], bp_pos);
            int tgt = (seg < 0) ? 0 : instructions[i].targets[seg];
            double tmr = (seg < 0) ? 0.0 : instructions[i].tmrcas[seg];
            if (same && (tgt != prev_target[i] || tmr != prev_tmrca[i]))
                same = false;
            prev_target[i] = tgt;
            prev_tmrca[i] = tmr;
            if (tgt >= 0)
                edges.push_back({i, tgt, tmr});
        }

        if (same) {
            for (int k = 0; k <= n; k++)
                afs[k] += prev_afs_contrib[k] * span;
            continue;
        }

        for (int i = 0; i < n; i++) { uf_parent[i] = i; uf_rank[i] = 0; uf_size[i] = 1; uf_height[i] = 0.0; }
        std::sort(edges.begin(), edges.end(),
                  [](const LightEdge& a, const LightEdge& b) { return a.height < b.height; });

        std::fill(prev_afs_contrib.begin(), prev_afs_contrib.end(), 0.0);
        for (const auto& e : edges) {
            int rc = find(e.child);
            int rp = find(e.parent);
            if (rc != rp) {
                int sc = uf_size[rc], sp = uf_size[rp];
                double hc = uf_height[rc], hp = uf_height[rp];
                if (hc < e.height) prev_afs_contrib[sc] += e.height - hc;
                if (hp < e.height) prev_afs_contrib[sp] += e.height - hp;
                int root = unite(rc, rp);
                uf_height[root] = e.height;
            }
        }
        for (int k = 0; k <= n; k++)
            afs[k] += prev_afs_contrib[k] * span;
    }
    return afs;
}


ThreadingInstructions::AssociationResult
ThreadingInstructions::association_diploid(const std::vector<double>& phenotypes,
                                            double p_threshold) const
{
    const int n = num_samples;
    const int n_dip = n / 2;
    AssociationResult result;
    if (n == 0 || num_sites == 0 || static_cast<int>(phenotypes.size()) != n_dip)
        return result;

    // Precompute phenotype stats
    double y_sum = 0.0, y_sum2 = 0.0;
    for (int i = 0; i < n_dip; i++) { y_sum += phenotypes[i]; y_sum2 += phenotypes[i] * phenotypes[i]; }
    double y_mean = y_sum / n_dip;
    double y_var = y_sum2 / n_dip - y_mean * y_mean;
    if (y_var <= 0.0) return result;

    // chi2 threshold from p-value (1 dof): chi2 = qchisq(1-p).
    // Use Abramowitz & Stegun approximation for the inverse normal, then square it.
    // For p = 5e-8, chi2 ~ 29.72
    double log_p = std::log(p_threshold);
    // Wilson-Hilferty approximation for chi2 inverse CDF (1 dof):
    // z = -sqrt(-2 * log(p/2))  (one-sided normal quantile)
    double z = -std::sqrt(-2.0 * std::log(p_threshold * 0.5));
    // Refine: Rational approximation for upper-tail normal quantile
    double t = std::sqrt(-2.0 * std::log(p_threshold * 0.5));
    z = t - (2.515517 + 0.802853 * t + 0.010328 * t * t) /
            (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t);
    double chi2_thresh = z * z;

    // Diploid phenotype sums indexed by haploid sample: pheno_val[hap] = phenotypes[hap/2]
    std::vector<double> pheno_dip(n_dip);
    for (int i = 0; i < n_dip; i++) pheno_dip[i] = phenotypes[i];

    // Use visit_clades to test each clade
    visit_clades([&](const std::vector<int>& descendants, double height, int start_bp, int end_bp) {
        int nd = static_cast<int>(descendants.size());
        if (nd < 2 || nd >= n - 1) return; // skip trivial clades

        // Sum phenotype over diploid individuals with at least one haploid in clade
        // For diploid: dosage = count of haploids in clade per individual
        // Simple approach: sum dosage * phenotype
        double dosage_sum = 0.0;
        double dosage_sum2 = 0.0;
        int dosage_count = 0;
        // Mark which diploid individuals have haploids in this clade
        // Use dosage: for each diploid i, dosage = #{hap in clade : hap/2 == i}
        std::vector<int> dosage(n_dip, 0);
        for (int h : descendants) dosage[h / 2]++;

        double gsum = 0.0, g2sum = 0.0, gy_sum = 0.0;
        for (int i = 0; i < n_dip; i++) {
            double g = dosage[i];
            gsum += g;
            g2sum += g * g;
            gy_sum += g * pheno_dip[i];
        }
        double g_mean = gsum / n_dip;
        double g_var = g2sum / n_dip - g_mean * g_mean;
        if (g_var <= 0.0) return;

        // Pearson correlation -> chi2 with 1 dof
        double cov = gy_sum / n_dip - g_mean * y_mean;
        double r2 = (cov * cov) / (g_var * y_var);
        double chi2 = r2 * n_dip;

        if (chi2 >= chi2_thresh) {
            // p-value: survival function of chi2(1) = 2 * Phi(-sqrt(chi2))
            // Use erfc for numerical stability
            double pval = std::erfc(std::sqrt(chi2 * 0.5));
            int midpoint = (start_bp + end_bp) / 2;
            result.position.push_back(midpoint);
            result.chi2.push_back(chi2);
            result.pvalue.push_back(pval);
            result.n_descendants.push_back(nd);
        }
    });

    return result;
}


ThreadingInstructions::MutationResult
ThreadingInstructions::generate_mutations(double mutation_rate, int seed) const
{
    MutationResult result;
    result.n_mutations = 0;
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return result;

    std::mt19937 rng(seed);

    visit_branches([&](const std::vector<int>& descendants, double child_h, double parent_h,
                        int start_bp, int end_bp) {
        double branch_length = parent_h - child_h;
        double span = static_cast<double>(end_bp - start_bp);
        double expected = mutation_rate * branch_length * span;
        std::poisson_distribution<int> pois(expected);
        int n_mut = pois(rng);

        std::uniform_real_distribution<double> pos_dist(start_bp, end_bp);
        for (int m = 0; m < n_mut; m++) {
            int pos = static_cast<int>(pos_dist(rng));
            result.positions.push_back(pos);
            // Build genotype row: 1 for descendants, 0 otherwise
            int offset = result.n_mutations * n;
            result.genotypes.resize(offset + n, 0);
            for (int d : descendants)
                result.genotypes[offset + d] = 1;
            result.n_mutations++;
        }
    });

    return result;
}


// ════════════════════════════════════════════════════════════════════════════
//  DAG (mutation-set) multiply — region-chunked for embarrassing parallelism
// ════════════════════════════════════════════════════════════════════════════

// Getters that aggregate across chunks
int ThreadingInstructions::get_dag_num_sets() const {
    int total = 0;
    for (const auto& c : dag_chunks) total += c.num_sets;
    return total;
}
int ThreadingInstructions::get_dag_total_edges() const {
    if (!dag_ready) return 0;
    int total = 0;
    for (const auto& c : dag_chunks) total += static_cast<int>(c.d2a.size());
    return total;
}
int ThreadingInstructions::get_dag_total_connections() const {
    if (!dag_ready) return 0;
    int total = 0;
    for (const auto& c : dag_chunks) total += static_cast<int>(c.indiv_set.size());
    return total;
}
int ThreadingInstructions::get_dag_total_mutations() const {
    if (!dag_ready) return 0;
    int total = 0;
    for (const auto& c : dag_chunks) total += static_cast<int>(c.mut_site.size());
    return total;
}

void ThreadingInstructions::prepare_dag_multiply(int num_chunks) {
    if (dag_ready) return;

    const int n = num_samples;
    const int m = num_sites;

    // Ensure ref genome is available for genotype_at_site()
    if (tree_ref_genome.empty()) {
        tree_ref_genome.assign(m, 0);
        for (int idx : instructions[0].mismatches) {
            if (idx >= 0 && idx < m) tree_ref_genome[idx] = 1;
        }
    }

    // Determine number of chunks
    if (num_chunks <= 0) {
        #ifdef _OPENMP
        num_chunks = omp_get_max_threads();
        #else
        num_chunks = 1;
        #endif
    }
    if (num_chunks > m) num_chunks = m;
    if (num_chunks < 1) num_chunks = 1;

    // ── Phase 1: Compute aligned split points (shared) ─────────────────────

    std::vector<std::set<int>> sample_splits(n);

    for (int i = 0; i < n; i++) {
        for (int s : instructions[i].starts) {
            sample_splits[i].insert(s);
        }
        sample_splits[i].insert(end);
        if (instructions[i].starts.empty() || instructions[i].starts[0] != start) {
            sample_splits[i].insert(start);
        }
    }

    // Reverse pass: propagate child boundaries to parents
    for (int i = n - 1; i >= 1; i--) {
        const auto& starts_i = instructions[i].starts;
        const auto& targets_i = instructions[i].targets;
        const int num_segs = static_cast<int>(starts_i.size());

        for (int k = 0; k < num_segs; k++) {
            const int target = targets_i[k];
            if (target == i) continue;  // self-ref
            const int seg_start_bp = starts_i[k];
            const int seg_end_bp = (k < num_segs - 1) ? starts_i[k + 1] : end;

            auto lo = sample_splits[i].lower_bound(seg_start_bp);
            auto hi = sample_splits[i].upper_bound(seg_end_bp);
            for (auto it = lo; it != hi; ++it) {
                sample_splits[target].insert(*it);
            }
        }
    }

    // ── Build chunks in parallel ──────────────────────────────────────────
    dag_chunks.resize(num_chunks);
    #pragma omp parallel for schedule(static)
    for (int c = 0; c < num_chunks; c++) {
        int site_lo = c * m / num_chunks;
        int site_hi = (c + 1) * m / num_chunks;
        dag_chunks[c] = build_dag_chunk(sample_splits, site_lo, site_hi);
    }

    dag_ready = true;
}

// Build a single DAG chunk for sites in [chunk_site_lo, chunk_site_hi).
ThreadingInstructions::DagChunk ThreadingInstructions::build_dag_chunk(
    const std::vector<std::set<int>>& sample_splits,
    int chunk_site_lo, int chunk_site_hi)
{
    const int n = num_samples;
    DagChunk chunk;
    chunk.site_lo = chunk_site_lo;
    chunk.site_hi = chunk_site_hi;

    // Helper lambda to flatten vector-of-vectors into flat + offset arrays
    auto flatten_int = [](const std::vector<std::vector<int>>& vv,
                          std::vector<int>& flat, std::vector<int>& off) {
        const int sz = static_cast<int>(vv.size());
        off.resize(sz + 1);
        off[0] = 0;
        for (int j = 0; j < sz; j++)
            off[j + 1] = off[j] + static_cast<int>(vv[j].size());
        flat.resize(off[sz]);
        for (int j = 0; j < sz; j++)
            std::copy(vv[j].begin(), vv[j].end(), flat.begin() + off[j]);
    };

    // ── Phase 2: Build DAG (forward pass) with site range filter ───────────

    std::vector<std::map<int, std::vector<int>>> sample_split_sets(n);
    for (int sp : sample_splits[0]) {
        sample_split_sets[0][sp] = {};
    }

    std::vector<std::vector<int>> tmp_d2a;
    std::vector<std::vector<int>> tmp_muts;
    std::vector<std::vector<int8_t>> tmp_signs;
    int mut_set_id = 0;

    for (int i = 1; i < n; i++) {
        const auto& starts_i = instructions[i].starts;
        const auto& targets_i = instructions[i].targets;
        const int num_segs = static_cast<int>(starts_i.size());

        std::vector<int> mm_i(instructions[i].mismatches.begin(),
                              instructions[i].mismatches.end());
        std::sort(mm_i.begin(), mm_i.end());

        auto& si_map = sample_split_sets[i];
        for (int sp : sample_splits[i]) {
            si_map[sp] = {};
        }

        for (int k = 0; k < num_segs; k++) {
            const int target = targets_i[k];
            const int seg_start_bp = starts_i[k];
            const int seg_end_bp = (k < num_segs - 1) ? starts_i[k + 1] : end;

            std::vector<int> my_splits;
            {
                auto lo = sample_splits[i].lower_bound(seg_start_bp);
                auto hi = sample_splits[i].lower_bound(seg_end_bp);
                for (auto it = lo; it != hi; ++it) my_splits.push_back(*it);
            }

            const auto& target_map = sample_split_sets[target];

            for (int r = 0; r < static_cast<int>(my_splits.size()); r++) {
                const int sub_start_bp = my_splits[r];
                const int sub_end_bp = (r + 1 < static_cast<int>(my_splits.size()))
                                       ? my_splits[r + 1] : seg_end_bp;
                if (sub_start_bp >= sub_end_bp) continue;

                std::vector<int> parent_sets;
                int count_meaningful = 0;
                {
                    auto t_lo = target_map.lower_bound(sub_start_bp);
                    auto t_hi = target_map.lower_bound(sub_end_bp);
                    for (auto t_it = t_lo; t_it != t_hi; ++t_it) {
                        if (!t_it->second.empty()) count_meaningful++;
                        parent_sets.insert(parent_sets.end(),
                                           t_it->second.begin(), t_it->second.end());
                    }
                }

                // Find mismatches intersected with this chunk's site range
                const int site_lo = static_cast<int>(
                    std::lower_bound(positions.begin(), positions.end(), sub_start_bp)
                    - positions.begin());
                const int site_hi = static_cast<int>(
                    std::lower_bound(positions.begin(), positions.end(), sub_end_bp)
                    - positions.begin());
                const int eff_lo = std::max(site_lo, chunk_site_lo);
                const int eff_hi = std::min(site_hi, chunk_site_hi);
                auto mm_lo = std::lower_bound(mm_i.begin(), mm_i.end(), eff_lo);
                auto mm_hi = std::lower_bound(mm_i.begin(), mm_i.end(), eff_hi);

                const bool has_mm = (eff_lo < eff_hi) && (mm_lo != mm_hi);

                if (!has_mm && count_meaningful <= 1) {
                    si_map[sub_start_bp] = std::move(parent_sets);
                } else {
                    std::vector<int> mm_in_range(mm_lo, mm_hi);
                    std::vector<int8_t> signs;
                    signs.reserve(mm_in_range.size());
                    for (int s : mm_in_range) {
                        const int g_target = genotype_at_site(target, s);
                        signs.push_back(static_cast<int8_t>(1 - 2 * g_target));
                    }

                    tmp_d2a.push_back(std::move(parent_sets));
                    tmp_muts.push_back(std::move(mm_in_range));
                    tmp_signs.push_back(std::move(signs));
                    si_map[sub_start_bp] = {mut_set_id};
                    mut_set_id++;
                }
            }
        }
    }

    chunk.num_sets = mut_set_id;

    if (chunk.num_sets == 0) {
        // Empty chunk: only ref_ones contribute
        chunk.d2a_off.assign(1, 0);
        chunk.a2d_off.assign(1, 0);
        chunk.mut_off.assign(1, 0);
        chunk.indiv_off.assign(n + 1, 0);
        chunk.set_indiv_off.assign(1, 0);
        for (int s = chunk_site_lo; s < chunk_site_hi; s++)
            if (tree_ref_genome[s]) chunk.ref_ones.push_back(s);
        return chunk;
    }

    // ── Phase 3: anc_to_desc + topological order ───────────────────────────

    std::vector<std::vector<int>> tmp_a2d(chunk.num_sets);
    for (int j = 0; j < chunk.num_sets; j++) {
        for (int anc : tmp_d2a[j]) {
            tmp_a2d[anc].push_back(j);
        }
    }

    std::vector<int> desc_count(chunk.num_sets);
    for (int j = 0; j < chunk.num_sets; j++)
        desc_count[j] = static_cast<int>(tmp_a2d[j].size());
    std::deque<int> queue;
    for (int j = 0; j < chunk.num_sets; j++)
        if (desc_count[j] == 0) queue.push_back(j);
    std::vector<int> topo_order;
    topo_order.reserve(chunk.num_sets);
    while (!queue.empty()) {
        int j = queue.front();
        queue.pop_front();
        topo_order.push_back(j);
        for (int anc : tmp_d2a[j])
            if (--desc_count[anc] == 0) queue.push_back(anc);
    }

    // ── Phase 4: Leaf connections ──────────────────────────────────────────

    std::vector<std::vector<int>> tmp_i2s(n);
    std::vector<std::vector<int>> tmp_s2i(chunk.num_sets);
    for (int i = 0; i < n; i++) {
        std::set<int> seen;
        for (const auto& kv : sample_split_sets[i]) {
            for (int s : kv.second) {
                if (seen.insert(s).second) {
                    tmp_i2s[i].push_back(s);
                    tmp_s2i[s].push_back(i);
                }
            }
        }
    }

    // ── Phase 5: Relabel set IDs to topo order ─────────────────────────────

    std::vector<int> relabel(chunk.num_sets);
    for (int ti = 0; ti < chunk.num_sets; ti++)
        relabel[topo_order[ti]] = ti;

    auto reorder = [&](std::vector<std::vector<int>>& vv, bool relabel_values) {
        std::vector<std::vector<int>> tmp(chunk.num_sets);
        for (int old_j = 0; old_j < chunk.num_sets; old_j++) {
            if (relabel_values)
                for (auto& v : vv[old_j]) v = relabel[v];
            tmp[relabel[old_j]] = std::move(vv[old_j]);
        }
        vv = std::move(tmp);
    };
    reorder(tmp_d2a, true);
    reorder(tmp_a2d, true);
    {
        std::vector<std::vector<int>> new_muts(chunk.num_sets);
        std::vector<std::vector<int8_t>> new_signs(chunk.num_sets);
        for (int old_j = 0; old_j < chunk.num_sets; old_j++) {
            new_muts[relabel[old_j]] = std::move(tmp_muts[old_j]);
            new_signs[relabel[old_j]] = std::move(tmp_signs[old_j]);
        }
        tmp_muts = std::move(new_muts);
        tmp_signs = std::move(new_signs);
    }
    for (int i = 0; i < n; i++)
        for (auto& s : tmp_i2s[i]) s = relabel[s];
    {
        std::vector<std::vector<int>> new_s2i(chunk.num_sets);
        for (int old_j = 0; old_j < chunk.num_sets; old_j++)
            new_s2i[relabel[old_j]] = std::move(tmp_s2i[old_j]);
        tmp_s2i = std::move(new_s2i);
    }

    // ── Phase 6: Flatten into offset arrays ────────────────────────────────

    flatten_int(tmp_d2a, chunk.d2a, chunk.d2a_off);
    flatten_int(tmp_a2d, chunk.a2d, chunk.a2d_off);
    flatten_int(tmp_i2s, chunk.indiv_set, chunk.indiv_off);
    flatten_int(tmp_s2i, chunk.set_indiv, chunk.set_indiv_off);

    chunk.mut_off.resize(chunk.num_sets + 1);
    chunk.mut_off[0] = 0;
    for (int j = 0; j < chunk.num_sets; j++)
        chunk.mut_off[j + 1] = chunk.mut_off[j] + static_cast<int>(tmp_muts[j].size());
    const int total_muts = chunk.mut_off[chunk.num_sets];
    chunk.mut_site.resize(total_muts);
    chunk.mut_sign.resize(total_muts);
    for (int j = 0; j < chunk.num_sets; j++) {
        std::copy(tmp_muts[j].begin(), tmp_muts[j].end(),
                  chunk.mut_site.begin() + chunk.mut_off[j]);
        std::copy(tmp_signs[j].begin(), tmp_signs[j].end(),
                  chunk.mut_sign.begin() + chunk.mut_off[j]);
    }

    for (int s = chunk_site_lo; s < chunk_site_hi; s++)
        if (tree_ref_genome[s]) chunk.ref_ones.push_back(s);

    return chunk;
}


// ── Right multiply single-column helper (no OMP — called from parallel regions) ──

void ThreadingInstructions::right_multiply_dag_single(const double* xp, double* out) const {
    const int n = num_samples;
    const int nc = static_cast<int>(dag_chunks.size());

    double total_ref_sum = 0.0;
    // Use a single flat partial buffer to avoid per-chunk allocations.
    // Chunk c's partials start at partial_base[c].
    std::vector<int> partial_base(nc + 1, 0);
    for (int c = 0; c < nc; c++) partial_base[c + 1] = partial_base[c] + dag_chunks[c].num_sets;
    std::vector<double> partial(partial_base[nc], 0.0);

    for (int c = 0; c < nc; c++) {
        const auto& ch = dag_chunks[c];
        double ref_sum = 0.0;
        for (int idx : ch.ref_ones) ref_sum += xp[idx];
        total_ref_sum += ref_sum;

        if (ch.num_sets > 0) {
            double* p = partial.data() + partial_base[c];
            for (int j = ch.num_sets - 1; j >= 0; j--) {
                double val = 0.0;
                for (int a = ch.d2a_off[j]; a < ch.d2a_off[j + 1]; a++)
                    val += p[ch.d2a[a]];
                for (int a = ch.mut_off[j]; a < ch.mut_off[j + 1]; a++)
                    val += ch.mut_sign[a] * xp[ch.mut_site[a]];
                p[j] = val;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        double sum = total_ref_sum;
        for (int c = 0; c < nc; c++) {
            const auto& ch = dag_chunks[c];
            if (ch.num_sets > 0) {
                const double* p = partial.data() + partial_base[c];
                for (int a = ch.indiv_off[i]; a < ch.indiv_off[i + 1]; a++)
                    sum += p[ch.indiv_set[a]];
            }
        }
        out[i] = sum;
    }
}


// ── Right multiply: result = G x  (single vector) ───────────────────────

std::vector<double> ThreadingInstructions::right_multiply_dag(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_sites)
        throw std::runtime_error("Input vector must have length " + std::to_string(num_sites));
    prepare_dag_multiply();

    std::vector<double> result(num_samples);
    right_multiply_dag_single(x.data(), result.data());
    return result;
}


// ── Right multiply batch: result = G X  (k vectors) ─────────────────────
// Hybrid parallelism: columns are split across threads in groups.
// Each thread processes its group with amortized multi-column DAG traversal
// (vectorization across columns) while different groups run in parallel.
// DAG structure is shared read-only in L3; each thread has private partials.

std::vector<double> ThreadingInstructions::right_multiply_dag_batch(
        const std::vector<double>& x_flat, int k) {
    if (static_cast<int>(x_flat.size()) != num_sites * k)
        throw std::runtime_error("Input must have length " + std::to_string(num_sites) + " * " + std::to_string(k));
    prepare_dag_multiply();

    const int n = num_samples;
    const int nc = static_cast<int>(dag_chunks.size());
    const double* xp = x_flat.data();

    std::vector<double> result(static_cast<size_t>(n) * k, 0.0);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int T = omp_get_num_threads();
        int cols_per = (k + T - 1) / T;
        int col_lo = std::min(tid * cols_per, k);
        int col_hi = std::min(col_lo + cols_per, k);
        int my_k = col_hi - col_lo;

        if (my_k > 0) {
            for (int c = 0; c < nc; c++) {
                const auto& ch = dag_chunks[c];

                // ref_sum for my columns
                std::vector<double> ref_sum(my_k, 0.0);
                for (int idx : ch.ref_ones) {
                    const double* xs = xp + static_cast<size_t>(idx) * k + col_lo;
                    for (int j = 0; j < my_k; j++) ref_sum[j] += xs[j];
                }

                if (ch.num_sets > 0) {
                    // Propagation (amortized across my_k columns)
                    std::vector<double> partial(static_cast<size_t>(ch.num_sets) * my_k, 0.0);
                    for (int j = ch.num_sets - 1; j >= 0; j--) {
                        double* pj = partial.data() + static_cast<size_t>(j) * my_k;
                        for (int a = ch.d2a_off[j]; a < ch.d2a_off[j + 1]; a++) {
                            const double* pa = partial.data() + static_cast<size_t>(ch.d2a[a]) * my_k;
                            for (int j2 = 0; j2 < my_k; j2++) pj[j2] += pa[j2];
                        }
                        for (int a = ch.mut_off[j]; a < ch.mut_off[j + 1]; a++) {
                            const int8_t sign = ch.mut_sign[a];
                            const double* xs = xp + static_cast<size_t>(ch.mut_site[a]) * k + col_lo;
                            for (int j2 = 0; j2 < my_k; j2++) pj[j2] += sign * xs[j2];
                        }
                    }

                    // Scatter — each thread writes to its own disjoint columns
                    for (int i = 0; i < n; i++) {
                        double* ri = result.data() + static_cast<size_t>(i) * k + col_lo;
                        for (int j = 0; j < my_k; j++) ri[j] += ref_sum[j];
                        for (int a = ch.indiv_off[i]; a < ch.indiv_off[i + 1]; a++) {
                            const double* pa = partial.data() + static_cast<size_t>(ch.indiv_set[a]) * my_k;
                            for (int j2 = 0; j2 < my_k; j2++) ri[j2] += pa[j2];
                        }
                    }
                } else {
                    for (int i = 0; i < n; i++) {
                        double* ri = result.data() + static_cast<size_t>(i) * k + col_lo;
                        for (int j = 0; j < my_k; j++) ri[j] += ref_sum[j];
                    }
                }
            }
        }
    }
    return result;
}


// ── Left multiply single-column helper (no OMP) ─────────────────────────

void ThreadingInstructions::left_multiply_dag_single(const double* yp, double* out) const {
    const int n = num_samples;
    const int nc = static_cast<int>(dag_chunks.size());

    double y_sum = 0.0;
    for (int i = 0; i < n; i++) y_sum += yp[i];

    // Zero output (caller may not have initialized)
    for (int s = 0; s < num_sites; s++) out[s] = 0.0;

    for (int c = 0; c < nc; c++) {
        const auto& ch = dag_chunks[c];

        for (int idx : ch.ref_ones) out[idx] = y_sum;

        if (ch.num_sets > 0) {
            std::vector<double> partial(ch.num_sets, 0.0);
            for (int j = 0; j < ch.num_sets; j++) {
                double sum = 0.0;
                for (int a = ch.set_indiv_off[j]; a < ch.set_indiv_off[j + 1]; a++)
                    sum += yp[ch.set_indiv[a]];
                partial[j] = sum;
            }
            for (int j = 0; j < ch.num_sets; j++) {
                for (int a = ch.a2d_off[j]; a < ch.a2d_off[j + 1]; a++)
                    partial[j] += partial[ch.a2d[a]];
            }
            for (int j = 0; j < ch.num_sets; j++) {
                const double pj = partial[j];
                for (int a = ch.mut_off[j]; a < ch.mut_off[j + 1]; a++)
                    out[ch.mut_site[a]] += ch.mut_sign[a] * pj;
            }
        }
    }
}


// ── Left multiply: result = G^T y  (single vector) ──────────────────────

std::vector<double> ThreadingInstructions::left_multiply_dag(const std::vector<double>& y) {
    if (static_cast<int>(y.size()) != num_samples)
        throw std::runtime_error("Input vector must have length " + std::to_string(num_samples));
    prepare_dag_multiply();

    std::vector<double> result(num_sites);
    left_multiply_dag_single(y.data(), result.data());
    return result;
}


// ── Left multiply batch: result = G^T Y  (k vectors) ────────────────────
// Hybrid parallelism: columns split across threads, amortized DAG traversal.

std::vector<double> ThreadingInstructions::left_multiply_dag_batch(
        const std::vector<double>& y_flat, int k) {
    if (static_cast<int>(y_flat.size()) != num_samples * k)
        throw std::runtime_error("Input must have length " + std::to_string(num_samples) + " * " + std::to_string(k));
    prepare_dag_multiply();

    const int n = num_samples;
    const int m = num_sites;
    const double* yp = y_flat.data();
    const int nc = static_cast<int>(dag_chunks.size());

    std::vector<double> result(static_cast<size_t>(m) * k, 0.0);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int T = omp_get_num_threads();
        int cols_per = (k + T - 1) / T;
        int col_lo = std::min(tid * cols_per, k);
        int col_hi = std::min(col_lo + cols_per, k);
        int my_k = col_hi - col_lo;

        if (my_k > 0) {
            // y_sum for my columns
            std::vector<double> y_sum(my_k, 0.0);
            for (int i = 0; i < n; i++) {
                const double* yi = yp + static_cast<size_t>(i) * k + col_lo;
                for (int j = 0; j < my_k; j++) y_sum[j] += yi[j];
            }

            for (int c = 0; c < nc; c++) {
                const auto& ch = dag_chunks[c];

                // ref sites get y_sum
                for (int idx : ch.ref_ones) {
                    double* rs = result.data() + static_cast<size_t>(idx) * k + col_lo;
                    for (int j = 0; j < my_k; j++) rs[j] = y_sum[j];
                }

                if (ch.num_sets > 0) {
                    // Gather from individuals (amortized across my_k columns)
                    std::vector<double> partial(static_cast<size_t>(ch.num_sets) * my_k, 0.0);
                    for (int j = 0; j < ch.num_sets; j++) {
                        double* pj = partial.data() + static_cast<size_t>(j) * my_k;
                        for (int a = ch.set_indiv_off[j]; a < ch.set_indiv_off[j + 1]; a++) {
                            const double* yi = yp + static_cast<size_t>(ch.set_indiv[a]) * k + col_lo;
                            for (int j2 = 0; j2 < my_k; j2++) pj[j2] += yi[j2];
                        }
                    }
                    // Propagation (forward/bottom-up)
                    for (int j = 0; j < ch.num_sets; j++) {
                        double* pj = partial.data() + static_cast<size_t>(j) * my_k;
                        for (int a = ch.a2d_off[j]; a < ch.a2d_off[j + 1]; a++) {
                            const double* pd = partial.data() + static_cast<size_t>(ch.a2d[a]) * my_k;
                            for (int j2 = 0; j2 < my_k; j2++) pj[j2] += pd[j2];
                        }
                    }
                    // Scatter to sites
                    for (int j = 0; j < ch.num_sets; j++) {
                        const double* pj = partial.data() + static_cast<size_t>(j) * my_k;
                        for (int a = ch.mut_off[j]; a < ch.mut_off[j + 1]; a++) {
                            const int s = ch.mut_site[a];
                            const int8_t sign = ch.mut_sign[a];
                            double* rs = result.data() + static_cast<size_t>(s) * k + col_lo;
                            for (int j2 = 0; j2 < my_k; j2++) rs[j2] += sign * pj[j2];
                        }
                    }
                }
            }
        }
    }
    return result;
}
