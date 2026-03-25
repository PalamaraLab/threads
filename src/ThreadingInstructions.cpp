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

    tree_ivl_seg.resize(static_cast<size_t>(n) * ni);
    for (int i = 0; i < n; i++) {
        const auto& starts_i = instructions[i].starts;
        const int n_segs = static_cast<int>(starts_i.size());
        int seg = 0;
        for (int j = 0; j < ni; j++) {
            const int site = tree_ivl_start[j];
            if (site == 0) {
                seg = 0;
            } else {
                const int pos = positions[site];
                while (seg + 1 < n_segs && starts_i[seg + 1] <= pos) {
                    seg++;
                }
            }
            tree_ivl_seg[static_cast<size_t>(i) * ni + j] = seg;
        }
    }

    // ── Phase 2: Per-sample segment-to-interval mapping ─────────────────
    // For each sample, record the first interval of each segment.
    // This enables segment-level iteration in the multiply.

    tree_seg_offset.resize(n);
    tree_seg_first_ivl.clear();
    for (int i = 0; i < n; i++) {
        tree_seg_offset[i] = static_cast<int>(tree_seg_first_ivl.size());
        int prev_seg = -1;
        for (int j = 0; j < ni; j++) {
            const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + j];
            if (seg != prev_seg) {
                tree_seg_first_ivl.push_back(j);
                prev_seg = seg;
            }
        }
    }
    // Sentinel for end-of-last-sample
    tree_seg_first_ivl.push_back(0);

    // ── Phase 3: Precompute mismatch signs and carry values ─────────────
    // Uses temporary reference-counted genotype cache.
    // Peak memory: O(m * max_active_targets) instead of O(n * m).

    // Count references: how many future samples use each target
    std::vector<int> ref_count(n, 0);
    for (int i = 1; i < n; i++) {
        std::set<int> seen;
        for (int t : instructions[i].targets) {
            if (t != i && seen.insert(t).second) ref_count[t]++;
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
    tree_mm_ivl.resize(total_mm);

    // Allocate carry storage
    tree_carry.assign(static_cast<size_t>(n) * ni, 0);

    // Reference-counted genotype cache
    std::unordered_map<int, std::vector<uint8_t>> geno_cache;
    std::vector<uint8_t> geno_buf(m, 0);

    for (int i = 0; i < n; i++) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());

        if (i == 0) {
            // Sample 0: genotype = ref_genome
            for (int s = 0; s < m; s++) geno_buf[s] = static_cast<uint8_t>(tree_ref_genome[s]);
        } else {
            // Build genotype from target + mismatches
            int carry = 0;
            for (int j = 0; j < ni; j++) {
                const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + j];
                const int target = instructions[i].targets[seg];
                const int a = tree_ivl_start[j];
                const int b = tree_ivl_end[j];

                if (target != i) {
                    const auto& tgt = geno_cache.at(target);
                    std::memcpy(&geno_buf[a], &tgt[a], b - a);
                    // Flip at mismatches
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, a);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, b);
                    for (const int* it = lo; it != hi; ++it) {
                        geno_buf[*it] ^= 1;
                    }
                    carry = geno_buf[b - 1];
                } else {
                    // Self-ref: carry forward, flip at mismatches
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, a);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, b);
                    const int* it = lo;
                    for (int s = a; s < b; s++) {
                        if (it != hi && *it == s) {
                            carry ^= 1;
                            ++it;
                        }
                        geno_buf[s] = static_cast<uint8_t>(carry);
                    }
                }
            }

            // Compute mismatch correction signs and interval mapping
            for (int k = 0; k < n_mm; k++) {
                const int s = mm_data[k];
                // Find interval containing site s (precompute for multiply)
                auto ivl_it = std::upper_bound(tree_ivl_start.begin(), tree_ivl_start.end(), s);
                const int j = static_cast<int>(ivl_it - tree_ivl_start.begin()) - 1;
                tree_mm_ivl[tree_mm_offset[i] + k] = j;
                const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + j];
                const int target = instructions[i].targets[seg];
                if (target != i) {
                    const int g_target = geno_cache.at(target)[s];
                    tree_mm_sign[tree_mm_offset[i] + k] = static_cast<int8_t>(1 - 2 * g_target);
                }
                // Self-ref mismatches: leave as 0 (not used)
            }
        }

        // Store carry at end of each interval
        for (int j = 0; j < ni; j++) {
            tree_carry[static_cast<size_t>(i) * ni + j] = static_cast<int8_t>(geno_buf[tree_ivl_end[j] - 1]);
        }

        // Cache genotype if still referenced
        if (ref_count[i] > 0) {
            geno_cache[i] = geno_buf;
        }

        // Decrement references and free unreferenced caches
        if (i > 0) {
            std::set<int> seen;
            for (int t : instructions[i].targets) {
                if (t != i && seen.insert(t).second) {
                    if (--ref_count[t] == 0) {
                        geno_cache.erase(t);
                    }
                }
            }
        }
    }

    tree_ready = true;
}

std::vector<double> ThreadingInstructions::right_multiply_tree(const std::vector<double>& x) {
    if (static_cast<int>(x.size()) != num_sites) {
        std::ostringstream oss;
        oss << "Input vector must have length " << num_sites << ".";
        throw std::runtime_error(oss.str());
    }

    prepare_tree_multiply();

    const int n = num_samples;
    const int m = num_sites;
    const int ni = tree_n_intervals;
    const double* xp = x.data();
    const int* ref = tree_ref_genome.data();
    const int* ivl_s = tree_ivl_start.data();
    const int* ivl_e = tree_ivl_end.data();

    // Prefix sums of x (for self-ref segments)
    std::vector<double> prefix_x(m + 1, 0.0);
    for (int s = 0; s < m; s++) {
        prefix_x[s + 1] = prefix_x[s] + xp[s];
    }

    // ── Reference-counted cumulative interval sums ──
    // Instead of n rows, allocate only O(tree_depth) rows via a pool.
    const size_t cum_stride = static_cast<size_t>(ni + 1);

    // Reference counting: how many future samples need each sample's cum row
    std::vector<int> cum_ref(n, 0);
    {
        std::vector<bool> seen(n, false);
        for (int i = 1; i < n; i++) {
            for (int t : instructions[i].targets)
                if (t != i && !seen[t]) { seen[t] = true; cum_ref[t]++; }
            for (int t : instructions[i].targets)
                if (t != i) seen[t] = false;
        }
    }

    // Pool of cum rows
    std::vector<std::vector<double>> cum_pool;
    std::vector<int> free_rows;
    std::vector<int> sample_to_row(n, -1);

    auto alloc_row = [&]() -> int {
        if (!free_rows.empty()) {
            int r = free_rows.back();
            free_rows.pop_back();
            std::fill(cum_pool[r].begin(), cum_pool[r].end(), 0.0);
            return r;
        }
        int r = static_cast<int>(cum_pool.size());
        cum_pool.emplace_back(cum_stride, 0.0);
        return r;
    };

    auto release_row = [&](int sample) {
        int r = sample_to_row[sample];
        if (r >= 0) {
            free_rows.push_back(r);
            sample_to_row[sample] = -1;
        }
    };

    // Sample 0: reference genome
    {
        int r0 = alloc_row();
        sample_to_row[0] = r0;
        double* c0 = cum_pool[r0].data();
        for (int j = 0; j < ni; j++) {
            double s = 0.0;
            for (int site = ivl_s[j]; site < ivl_e[j]; site++) {
                s += xp[site] * ref[site];
            }
            c0[j + 1] = c0[j] + s;
        }
    }

    std::vector<double> out(n, 0.0);
    out[0] = cum_pool[sample_to_row[0]][ni];

    // Process samples 1..n-1 using segment-level loop
    for (int i = 1; i < n; i++) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        const bool need_cum = (cum_ref[i] > 0);

        int my_row = -1;
        double* my_cum = nullptr;
        if (need_cum) {
            my_row = alloc_row();
            sample_to_row[i] = my_row;
            my_cum = cum_pool[my_row].data();
        }

        const int seg_off = tree_seg_offset[i];
        int n_mapped_segs;
        {
            int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                        : static_cast<int>(tree_seg_first_ivl.size()) - 1;
            n_mapped_segs = next_off - seg_off;
        }

        double total = 0.0;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int first_ivl = tree_seg_first_ivl[seg_off + sk];
            const int last_ivl = (sk + 1 < n_mapped_segs)
                                     ? tree_seg_first_ivl[seg_off + sk + 1]
                                     : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + first_ivl];
            const int target = instructions[i].targets[seg];
            const int site_a = ivl_s[first_ivl];
            const int site_b = ivl_e[last_ivl - 1];

            if (target != i) {
                const double* tgt_cum = cum_pool[sample_to_row[target]].data();
                const double base = tgt_cum[last_ivl] - tgt_cum[first_ivl];

                const int* lo = std::lower_bound(mm_data, mm_data + n_mm, site_a);
                const int* hi = std::lower_bound(mm_data, mm_data + n_mm, site_b);
                double corr = 0.0;
                for (const int* it = lo; it != hi; ++it) {
                    const int k = static_cast<int>(it - mm_data);
                    corr += xp[*it] * my_signs[k];
                }
                total += base + corr;

                if (need_cum) {
                    const double offset = my_cum[first_ivl] - tgt_cum[first_ivl];
                    std::memcpy(&my_cum[first_ivl + 1], &tgt_cum[first_ivl + 1],
                                (last_ivl - first_ivl) * sizeof(double));
                    for (int j = first_ivl + 1; j <= last_ivl; j++) {
                        my_cum[j] += offset;
                    }
                    const int* mm_ivl_data = tree_mm_ivl.data() + tree_mm_offset[i];
                    for (const int* it = lo; it != hi; ++it) {
                        const int k = static_cast<int>(it - mm_data);
                        const int mm_ivl = mm_ivl_data[k];
                        const double c = xp[*it] * my_signs[k];
                        for (int j = mm_ivl + 1; j <= last_ivl; j++) {
                            my_cum[j] += c;
                        }
                    }
                }
            } else {
                // Self-ref: interval by interval
                for (int j = first_ivl; j < last_ivl; j++) {
                    const int a = ivl_s[j];
                    const int b = ivl_e[j];
                    int carry = (j == 0) ? 0
                                : static_cast<int>(tree_carry[static_cast<size_t>(i) * ni + j - 1]);
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, a);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, b);
                    double ivl_sum = 0.0;
                    int prev = a;
                    for (const int* it = lo; it != hi; ++it) {
                        const int mm_site = *it;
                        if (mm_site > prev)
                            ivl_sum += carry * (prefix_x[mm_site] - prefix_x[prev]);
                        carry = 1 - carry;
                        ivl_sum += xp[mm_site] * carry;
                        prev = mm_site + 1;
                    }
                    if (b > prev)
                        ivl_sum += carry * (prefix_x[b] - prefix_x[prev]);
                    if (need_cum) my_cum[j + 1] = my_cum[j] + ivl_sum;
                    total += ivl_sum;
                }
            }
        }
        out[i] = total;

        // Release targets no longer needed
        {
            const auto& tgts = instructions[i].targets;
            const int nt = static_cast<int>(tgts.size());
            for (int si = 0; si < nt; si++) {
                int t = tgts[si];
                if (t == i) continue;
                bool dup = false;
                for (int sj = 0; sj < si; sj++)
                    if (tgts[sj] == t) { dup = true; break; }
                if (!dup && --cum_ref[t] == 0) release_row(t);
            }
        }
        // Release own row if no one references us
        if (cum_ref[i] == 0) release_row(i);
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

    // Per-sample accumulated weight matrix: W[i * ni + j] = weight for sample i at interval j
    // Initialized to x[i] for all intervals, then accumulated bottom-up.
    std::vector<double> W(static_cast<size_t>(n) * ni);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; i++) {
        const double xi = xp[i];
        double* row = W.data() + static_cast<size_t>(i) * ni;
        std::fill(row, row + ni, xi);
    }

    // ── Pass 1 (sequential): Bottom-up W propagation ───────────────────────
    // Push each sample's weight to its target. This is O(n * avg_intervals)
    // and must be sequential due to tree dependencies (child before parent).
    for (int i = n - 1; i >= 1; i--) {
        const double* Wi = &W[static_cast<size_t>(i) * ni];
        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped_segs = next_off - seg_off;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int first_ivl = tree_seg_first_ivl[seg_off + sk];
            const int last_ivl = (sk + 1 < n_mapped_segs)
                                     ? tree_seg_first_ivl[seg_off + sk + 1]
                                     : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + first_ivl];
            const int target = instructions[i].targets[seg];

            if (target != i) {
                double* Wt = &W[static_cast<size_t>(target) * ni];
                for (int j = first_ivl; j < last_ivl; j++) {
                    Wt[j] += Wi[j];
                }
            }
        }
    }

    // ── Pass 2 (parallel): Accumulate corrections into out[] ─────────────
    // W is now finalized. Each sample's mismatch corrections and self-ref
    // carry fills are independent. Use thread-local buffers to avoid contention.
    std::vector<double> out(m, 0.0);

#ifdef _OPENMP
    #pragma omp parallel
    {
        std::vector<double> local_out(m, 0.0);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 1; i < n; i++) {
            const int* mm_data = instructions[i].mismatches.data();
            const int n_mm = static_cast<int>(instructions[i].mismatches.size());
            const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
            const double* Wi = &W[static_cast<size_t>(i) * ni];

            const int seg_off = tree_seg_offset[i];
            const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                              : static_cast<int>(tree_seg_first_ivl.size()) - 1;
            const int n_mapped_segs = next_off - seg_off;

            for (int sk = 0; sk < n_mapped_segs; sk++) {
                const int first_ivl = tree_seg_first_ivl[seg_off + sk];
                const int last_ivl = (sk + 1 < n_mapped_segs)
                                         ? tree_seg_first_ivl[seg_off + sk + 1]
                                         : ni;
                const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + first_ivl];
                const int target = instructions[i].targets[seg];
                const int site_a = tree_ivl_start[first_ivl];
                const int site_b = tree_ivl_end[last_ivl - 1];

                if (target != i) {
                    // Mismatch corrections
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, site_a);
                    const int* hi = std::lower_bound(mm_data, mm_data + n_mm, site_b);
                    const int* mm_ivl_data = tree_mm_ivl.data() + tree_mm_offset[i];
                    for (const int* it = lo; it != hi; ++it) {
                        const int k = static_cast<int>(it - mm_data);
                        local_out[*it] += Wi[mm_ivl_data[k]] * my_signs[k];
                    }
                } else {
                    // Self-ref: carry-run decomposition
                    for (int j = first_ivl; j < last_ivl; j++) {
                        const int a = tree_ivl_start[j];
                        const int b = tree_ivl_end[j];
                        int carry = (j == 0) ? 0
                                    : static_cast<int>(tree_carry[static_cast<size_t>(i) * ni + j - 1]);
                        const int* lo = std::lower_bound(mm_data, mm_data + n_mm, a);
                        const int* hi_mm = std::lower_bound(mm_data, mm_data + n_mm, b);
                        const double wi = Wi[j];
                        int prev = a;
                        double* lop = local_out.data();
                        for (const int* it = lo; it != hi_mm; ++it) {
                            const int mm_site = *it;
                            if (carry) {
                                double* __restrict__ op = lop + prev;
                                const int len = mm_site - prev;
                                for (int s = 0; s < len; s++) op[s] += wi;
                            }
                            carry ^= 1;
                            lop[mm_site] += wi * carry;
                            prev = mm_site + 1;
                        }
                        if (carry && b > prev) {
                            double* __restrict__ op = lop + prev;
                            const int len = b - prev;
                            for (int s = 0; s < len; s++) op[s] += wi;
                        }
                    }
                }
            }
        }

        // Reduce thread-local buffers into shared out
        #pragma omp critical
        {
            double* outp = out.data();
            const double* lp = local_out.data();
            for (int s = 0; s < m; s++) outp[s] += lp[s];
        }
    }
#else
    // Single-threaded fallback
    for (int i = n - 1; i >= 1; i--) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        const double* Wi = &W[static_cast<size_t>(i) * ni];

        const int seg_off = tree_seg_offset[i];
        const int next_off = (i + 1 < n) ? tree_seg_offset[i + 1]
                                          : static_cast<int>(tree_seg_first_ivl.size()) - 1;
        const int n_mapped_segs = next_off - seg_off;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int first_ivl = tree_seg_first_ivl[seg_off + sk];
            const int last_ivl = (sk + 1 < n_mapped_segs)
                                     ? tree_seg_first_ivl[seg_off + sk + 1]
                                     : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i) * ni + first_ivl];
            const int target = instructions[i].targets[seg];
            const int site_a = tree_ivl_start[first_ivl];
            const int site_b = tree_ivl_end[last_ivl - 1];

            if (target != i) {
                const int* lo = std::lower_bound(mm_data, mm_data + n_mm, site_a);
                const int* hi = std::lower_bound(mm_data, mm_data + n_mm, site_b);
                const int* mm_ivl_data = tree_mm_ivl.data() + tree_mm_offset[i];
                for (const int* it = lo; it != hi; ++it) {
                    const int k = static_cast<int>(it - mm_data);
                    out[*it] += Wi[mm_ivl_data[k]] * my_signs[k];
                }
            } else {
                for (int j = first_ivl; j < last_ivl; j++) {
                    const int a = tree_ivl_start[j];
                    const int b = tree_ivl_end[j];
                    int carry = (j == 0) ? 0
                                : static_cast<int>(tree_carry[static_cast<size_t>(i) * ni + j - 1]);
                    const int* lo = std::lower_bound(mm_data, mm_data + n_mm, a);
                    const int* hi_mm = std::lower_bound(mm_data, mm_data + n_mm, b);
                    const double wi = Wi[j];
                    int prev = a;
                    double* outp = out.data();
                    for (const int* it = lo; it != hi_mm; ++it) {
                        const int mm_site = *it;
                        if (carry) {
                            double* __restrict__ op = outp + prev;
                            const int len = mm_site - prev;
                            for (int s = 0; s < len; s++) op[s] += wi;
                        }
                        carry ^= 1;
                        outp[mm_site] += wi * carry;
                        prev = mm_site + 1;
                    }
                    if (carry && b > prev) {
                        double* __restrict__ op = outp + prev;
                        const int len = b - prev;
                        for (int s = 0; s < len; s++) op[s] += wi;
                    }
                }
            }
        }
    }
#endif

    // Sample 0: multiply accumulated weight by reference genome
    {
        const double* W0 = &W[0];
        const int* ref = tree_ref_genome.data();
        for (int j = 0; j < ni; j++) {
            const double w0j = W0[j];
            for (int s = tree_ivl_start[j]; s < tree_ivl_end[j]; s++) {
                out[s] += w0j * ref[s];
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

    prepare_tree_multiply();
    const int n = num_samples, m = num_sites, ni = tree_n_intervals;
    const double* xp = x_flat.data();
    const int* ref = tree_ref_genome.data();
    const int* ivl_s = tree_ivl_start.data();
    const int* ivl_e = tree_ivl_end.data();

    // k prefix sums
    std::vector<double> prefix_x(static_cast<size_t>(m + 1) * k, 0.0);
    for (int s = 0; s < m; s++) {
        const double* xs = &xp[s * k];
        const double* ps = &prefix_x[s * k];
        double* pd = &prefix_x[(s + 1) * k];
        for (int c = 0; c < k; c++) pd[c] = ps[c] + xs[c];
    }

    const size_t cum_stride = static_cast<size_t>(ni + 1) * k;
    std::vector<int> cum_ref(n, 0);
    {
        std::vector<bool> seen(n, false);
        for (int i = 1; i < n; i++) {
            for (int t : instructions[i].targets)
                if (t != i && !seen[t]) { seen[t] = true; cum_ref[t]++; }
            for (int t : instructions[i].targets)
                if (t != i) seen[t] = false;
        }
    }
    std::vector<std::vector<double>> cum_pool;
    std::vector<int> free_rows, sample_to_row(n, -1);
    auto alloc_row = [&]() -> int {
        if (!free_rows.empty()) {
            int r = free_rows.back(); free_rows.pop_back();
            std::fill(cum_pool[r].begin(), cum_pool[r].end(), 0.0);
            return r;
        }
        int r = static_cast<int>(cum_pool.size());
        cum_pool.emplace_back(cum_stride, 0.0);
        return r;
    };
    auto release_row = [&](int s) {
        int r = sample_to_row[s];
        if (r >= 0) { free_rows.push_back(r); sample_to_row[s] = -1; }
    };

    // Sample 0
    { int r0 = alloc_row(); sample_to_row[0] = r0;
      double* c0 = cum_pool[r0].data();
      for (int j = 0; j < ni; j++) {
          double* dst = &c0[(j + 1) * k];
          const double* src = &c0[j * k];
          for (int c = 0; c < k; c++) dst[c] = src[c];
          for (int site = ivl_s[j]; site < ivl_e[j]; site++)
              if (ref[site]) { const double* xs = &xp[site*k];
                  for (int c = 0; c < k; c++) dst[c] += xs[c]; }
      }
    }

    std::vector<double> out(static_cast<size_t>(n) * k, 0.0);
    { const double* c0e = &cum_pool[sample_to_row[0]][ni * k];
      for (int c = 0; c < k; c++) out[c] = c0e[c]; }

    for (int i = 1; i < n; i++) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        const bool need_cum = (cum_ref[i] > 0);
        int my_row = -1;
        double* my_cum = nullptr;
        if (need_cum) {
            my_row = alloc_row(); sample_to_row[i] = my_row;
            my_cum = cum_pool[my_row].data();
        }
        const int seg_off = tree_seg_offset[i];
        int n_mapped_segs = ((i+1<n) ? tree_seg_offset[i+1]
                              : static_cast<int>(tree_seg_first_ivl.size())-1) - seg_off;
        double* sample_out = &out[i * k];

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int fi = tree_seg_first_ivl[seg_off + sk];
            const int li = (sk+1<n_mapped_segs) ? tree_seg_first_ivl[seg_off+sk+1] : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i)*ni + fi];
            const int target = instructions[i].targets[seg];
            const int sa = ivl_s[fi], sb = ivl_e[li-1];

            if (target != i) {
                const double* tc = cum_pool[sample_to_row[target]].data();
                const int* lo = std::lower_bound(mm_data, mm_data+n_mm, sa);
                const int* hi = std::lower_bound(mm_data, mm_data+n_mm, sb);

                // Compute output: base from target cum + mismatch correction
                for (int c = 0; c < k; c++) {
                    double base = tc[li*k+c] - tc[fi*k+c], corr = 0.0;
                    for (const int* it = lo; it != hi; ++it)
                        corr += xp[*it*k+c] * my_signs[it-mm_data];
                    sample_out[c] += base + corr;
                }
                // Build cum only if someone will reference us
                if (need_cum) {
                    const double* fi_cum = &my_cum[fi*k];
                    const double* fi_tc = &tc[fi*k];
                    for (int j = fi+1; j <= li; j++)
                        for (int c = 0; c < k; c++)
                            my_cum[j*k+c] = tc[j*k+c] + fi_cum[c] - fi_tc[c];
                    // Apply mismatch corrections via deferred deltas + prefix sum
                    const int n_seg_ivls = li - fi;
                    if (lo != hi && n_seg_ivls > 0) {
                        std::vector<double> mm_corr(static_cast<size_t>(n_seg_ivls) * k, 0.0);
                        const int* miv = tree_mm_ivl.data() + tree_mm_offset[i];
                        for (const int* it = lo; it != hi; ++it) {
                            int kk = static_cast<int>(it-mm_data);
                            int mivl = miv[kk]; double sv = my_signs[kk];
                            const double* xs = &xp[*it*k];
                            int idx = mivl - fi;
                            if (idx >= 0 && idx < n_seg_ivls) {
                                double* d = &mm_corr[idx * k];
                                for (int c = 0; c < k; c++) d[c] += xs[c]*sv;
                            }
                        }
                        double* running = &mm_corr[0];
                        for (int c = 0; c < k; c++)
                            my_cum[(fi+1)*k+c] += running[c];
                        for (int j = 1; j < n_seg_ivls; j++) {
                            double* cur = &mm_corr[j*k];
                            const double* prv = &mm_corr[(j-1)*k];
                            for (int c = 0; c < k; c++) {
                                cur[c] += prv[c];
                                my_cum[(fi+1+j)*k+c] += cur[c];
                            }
                        }
                    }
                }
            } else {
                for (int j = fi; j < li; j++) {
                    int a = ivl_s[j], b = ivl_e[j];
                    int carry = (j==0)?0:static_cast<int>(tree_carry[static_cast<size_t>(i)*ni+j-1]);
                    const int* lo = std::lower_bound(mm_data, mm_data+n_mm, a);
                    const int* hi = std::lower_bound(mm_data, mm_data+n_mm, b);
                    for (int c = 0; c < k; c++) {
                        double ivl_sum = 0.0; int cc = carry; int prev = a;
                        for (const int* it = lo; it != hi; ++it) {
                            int ms = *it;
                            if (ms > prev) ivl_sum += cc*(prefix_x[ms*k+c]-prefix_x[prev*k+c]);
                            cc = 1-cc; ivl_sum += xp[ms*k+c]*cc; prev = ms+1;
                        }
                        if (b > prev) ivl_sum += cc*(prefix_x[b*k+c]-prefix_x[prev*k+c]);
                        if (need_cum) my_cum[(j+1)*k+c] = my_cum[j*k+c]+ivl_sum;
                        sample_out[c] += ivl_sum;
                    }
                }
            }
        }
        {
            const auto& tgts = instructions[i].targets;
            const int nt = static_cast<int>(tgts.size());
            for (int si = 0; si < nt; si++) {
                int t = tgts[si];
                if (t == i) continue;
                bool dup = false;
                for (int sj = 0; sj < si; sj++)
                    if (tgts[sj] == t) { dup = true; break; }
                if (!dup && --cum_ref[t] == 0) release_row(t);
            }
        }
        if (cum_ref[i]==0) release_row(i);
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

    const size_t Wstride = static_cast<size_t>(ni) * k;
    std::vector<double> W(static_cast<size_t>(n) * Wstride);
    for (int i = 0; i < n; i++) {
        const double* xi = &xp[i * k];
        double* Wi = &W[i * Wstride];
        for (int j = 0; j < ni; j++) {
            double* d = &Wi[j*k];
            for (int c = 0; c < k; c++) d[c] = xi[c];
        }
    }

    std::vector<double> out(static_cast<size_t>(m) * k, 0.0);

    for (int i = n-1; i >= 1; i--) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        const double* Wi = &W[i * Wstride];
        const int seg_off = tree_seg_offset[i];
        int n_mapped_segs = ((i+1<n) ? tree_seg_offset[i+1]
                              : static_cast<int>(tree_seg_first_ivl.size())-1) - seg_off;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int fi = tree_seg_first_ivl[seg_off+sk];
            const int li = (sk+1<n_mapped_segs) ? tree_seg_first_ivl[seg_off+sk+1] : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i)*ni+fi];
            const int target = instructions[i].targets[seg];
            const int sa = tree_ivl_start[fi], sb = tree_ivl_end[li-1];

            if (target != i) {
                double* Wt = &W[static_cast<size_t>(target)*Wstride];
                for (int j = fi; j < li; j++) {
                    const double* s = &Wi[j*k]; double* d = &Wt[j*k];
                    for (int c = 0; c < k; c++) d[c] += s[c];
                }
                const int* lo = std::lower_bound(mm_data, mm_data+n_mm, sa);
                const int* hi = std::lower_bound(mm_data, mm_data+n_mm, sb);
                const int* miv = tree_mm_ivl.data() + tree_mm_offset[i];
                for (const int* it = lo; it != hi; ++it) {
                    int kk = static_cast<int>(it-mm_data);
                    double sv = my_signs[kk];
                    const double* wi = &Wi[miv[kk]*k];
                    double* od = &out[*it*k];
                    for (int c = 0; c < k; c++) od[c] += wi[c]*sv;
                }
            } else {
                for (int j = fi; j < li; j++) {
                    int a = tree_ivl_start[j], b = tree_ivl_end[j];
                    int carry = (j==0)?0:static_cast<int>(tree_carry[static_cast<size_t>(i)*ni+j-1]);
                    const int* lo = std::lower_bound(mm_data, mm_data+n_mm, a);
                    const int* hi_mm = std::lower_bound(mm_data, mm_data+n_mm, b);
                    const double* wi = &Wi[j*k];
                    int prev = a;
                    for (const int* it = lo; it != hi_mm; ++it) {
                        int ms = *it;
                        if (carry && ms > prev)
                            for (int s = prev; s < ms; s++) {
                                double* od = &out[s*k];
                                for (int c = 0; c < k; c++) od[c] += wi[c];
                            }
                        carry ^= 1;
                        if (carry) { double* od = &out[ms*k];
                            for (int c = 0; c < k; c++) od[c] += wi[c]; }
                        prev = ms+1;
                    }
                    if (carry && b > prev)
                        for (int s = prev; s < b; s++) {
                            double* od = &out[s*k];
                            for (int c = 0; c < k; c++) od[c] += wi[c];
                        }
                }
            }
        }
    }

    // Sample 0
    { const double* W0 = &W[0];
      const int* ref = tree_ref_genome.data();
      for (int j = 0; j < ni; j++) {
          const double* w0j = &W0[j*k];
          for (int s = tree_ivl_start[j]; s < tree_ivl_end[j]; s++)
              if (ref[s]) { double* od = &out[s*k];
                  for (int c = 0; c < k; c++) od[c] += w0j[c]; }
      }
    }
    return out;
}


// ═══════════════════════════════════════════════════════════════════════════
// Range-restricted tree multiply.
// Uses the FULL precomputed tree structure (correct carry state) but only
// accumulates contributions from sites in [site_start, site_end).
// ═══════════════════════════════════════════════════════════════════════════

std::vector<double> ThreadingInstructions::right_multiply_tree_range(
        const std::vector<double>& x, int site_start, int site_end) {
    const int range_len = site_end - site_start;
    if (static_cast<int>(x.size()) != range_len)
        throw std::runtime_error("Input length must equal site_end - site_start");

    prepare_tree_multiply();
    const int n = num_samples, m = num_sites, ni = tree_n_intervals;
    const int* ref = tree_ref_genome.data();
    const int* ivl_s = tree_ivl_start.data();
    const int* ivl_e = tree_ivl_end.data();

    // Build full-length x_full with zeros outside range
    std::vector<double> x_full(m, 0.0);
    for (int s = 0; s < range_len; s++) x_full[site_start + s] = x[s];
    const double* xp = x_full.data();

    // Prefix sums of x_full
    std::vector<double> prefix_x(m + 1, 0.0);
    for (int s = 0; s < m; s++) prefix_x[s + 1] = prefix_x[s] + xp[s];

    // Same cum-pool logic as right_multiply_tree
    const size_t cum_stride = static_cast<size_t>(ni + 1);
    std::vector<int> cum_ref(n, 0);
    for (int i = 1; i < n; i++) {
        std::set<int> seen;
        for (int t : instructions[i].targets)
            if (t != i && seen.insert(t).second) cum_ref[t]++;
    }
    std::vector<std::vector<double>> cum_pool;
    std::vector<int> free_rows, sample_to_row(n, -1);
    auto alloc_row = [&]() -> int {
        if (!free_rows.empty()) {
            int r = free_rows.back(); free_rows.pop_back();
            std::fill(cum_pool[r].begin(), cum_pool[r].end(), 0.0);
            return r;
        }
        int r = static_cast<int>(cum_pool.size());
        cum_pool.emplace_back(cum_stride, 0.0);
        return r;
    };
    auto release_row = [&](int s) {
        int r = sample_to_row[s];
        if (r >= 0) { free_rows.push_back(r); sample_to_row[s] = -1; }
    };

    { int r0 = alloc_row(); sample_to_row[0] = r0;
      double* c0 = cum_pool[r0].data();
      for (int j = 0; j < ni; j++) {
          double s = 0.0;
          for (int site = ivl_s[j]; site < ivl_e[j]; site++) s += xp[site] * ref[site];
          c0[j + 1] = c0[j] + s;
      }
    }

    std::vector<double> out(n, 0.0);
    out[0] = cum_pool[sample_to_row[0]][ni];

    for (int i = 1; i < n; i++) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        int my_row = alloc_row(); sample_to_row[i] = my_row;
        double* my_cum = cum_pool[my_row].data();
        const int seg_off = tree_seg_offset[i];
        int n_mapped_segs = ((i+1<n) ? tree_seg_offset[i+1]
                              : static_cast<int>(tree_seg_first_ivl.size())-1) - seg_off;
        double total = 0.0;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int fi = tree_seg_first_ivl[seg_off + sk];
            const int li = (sk+1<n_mapped_segs) ? tree_seg_first_ivl[seg_off+sk+1] : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i)*ni + fi];
            const int target = instructions[i].targets[seg];
            const int sa = ivl_s[fi], sb = ivl_e[li-1];

            if (target != i) {
                const double* tc = cum_pool[sample_to_row[target]].data();
                const double base = tc[li] - tc[fi];
                const int* lo = std::lower_bound(mm_data, mm_data+n_mm, sa);
                const int* hi = std::lower_bound(mm_data, mm_data+n_mm, sb);
                double corr = 0.0;
                for (const int* it = lo; it != hi; ++it)
                    corr += xp[*it] * my_signs[it-mm_data];
                total += base + corr;
                const double offset = my_cum[fi] - tc[fi];
                std::memcpy(&my_cum[fi+1], &tc[fi+1], (li-fi)*sizeof(double));
                for (int j = fi+1; j <= li; j++) my_cum[j] += offset;
                const int* miv = tree_mm_ivl.data() + tree_mm_offset[i];
                for (const int* it = lo; it != hi; ++it) {
                    int kk = static_cast<int>(it-mm_data);
                    double c = xp[*it] * my_signs[kk];
                    for (int j = miv[kk]+1; j <= li; j++) my_cum[j] += c;
                }
            } else {
                for (int j = fi; j < li; j++) {
                    int a = ivl_s[j], b = ivl_e[j];
                    int carry = (j==0)?0:static_cast<int>(tree_carry[static_cast<size_t>(i)*ni+j-1]);
                    const int* lo = std::lower_bound(mm_data, mm_data+n_mm, a);
                    const int* hi = std::lower_bound(mm_data, mm_data+n_mm, b);
                    double ivl_sum = 0.0; int prev = a;
                    for (const int* it = lo; it != hi; ++it) {
                        int ms = *it;
                        if (ms > prev) ivl_sum += carry*(prefix_x[ms]-prefix_x[prev]);
                        carry = 1-carry; ivl_sum += xp[ms]*carry; prev = ms+1;
                    }
                    if (b > prev) ivl_sum += carry*(prefix_x[b]-prefix_x[prev]);
                    my_cum[j+1] = my_cum[j]+ivl_sum; total += ivl_sum;
                }
            }
        }
        out[i] = total;
        { std::set<int> seen;
          for (int t : instructions[i].targets)
              if (t!=i && seen.insert(t).second && --cum_ref[t]==0) release_row(t);
        }
        if (cum_ref[i]==0) release_row(i);
    }
    return out;
}


std::vector<double> ThreadingInstructions::left_multiply_tree_range(
        const std::vector<double>& x, int site_start, int site_end) {
    if (static_cast<int>(x.size()) != num_samples)
        throw std::runtime_error("Input must have length num_samples");

    prepare_tree_multiply();
    const int n = num_samples, ni = tree_n_intervals;
    const double* xp = x.data();

    // Full W matrix needed for correct weight propagation
    std::vector<double> W(static_cast<size_t>(n) * ni);
    for (int i = 0; i < n; i++) {
        double* row = W.data() + static_cast<size_t>(i) * ni;
        std::fill(row, row + ni, xp[i]);
    }

    const int range_len = site_end - site_start;
    if (site_start < 0 || site_end > num_sites || site_start >= site_end)
        throw std::runtime_error("Invalid site range");
    std::vector<double> out(range_len, 0.0);

    // Bottom-up: full W propagation, restricted output
    for (int i = n-1; i >= 1; i--) {
        const int* mm_data = instructions[i].mismatches.data();
        const int n_mm = static_cast<int>(instructions[i].mismatches.size());
        const int8_t* my_signs = &tree_mm_sign[tree_mm_offset[i]];
        const double* Wi = &W[static_cast<size_t>(i) * ni];
        const int seg_off = tree_seg_offset[i];
        int n_mapped_segs = ((i+1<n) ? tree_seg_offset[i+1]
                              : static_cast<int>(tree_seg_first_ivl.size())-1) - seg_off;

        for (int sk = 0; sk < n_mapped_segs; sk++) {
            const int fi = tree_seg_first_ivl[seg_off+sk];
            const int li = (sk+1<n_mapped_segs) ? tree_seg_first_ivl[seg_off+sk+1] : ni;
            const int seg = tree_ivl_seg[static_cast<size_t>(i)*ni+fi];
            const int target = instructions[i].targets[seg];
            const int sa = tree_ivl_start[fi], sb = tree_ivl_end[li-1];

            if (target != i) {
                // W push: always needed for correct propagation
                double* Wt = &W[static_cast<size_t>(target) * ni];
                for (int j = fi; j < li; j++) Wt[j] += Wi[j];
                // Mismatch corrections: only for sites overlapping range
                if (sa < site_end && sb > site_start) {
                    const int range_lo = std::max(sa, site_start);
                    const int range_hi = std::min(sb, site_end);
                    const int* lo = std::lower_bound(mm_data, mm_data+n_mm, range_lo);
                    const int* hi = std::lower_bound(mm_data, mm_data+n_mm, range_hi);
                    const int* miv = tree_mm_ivl.data() + tree_mm_offset[i];
                    for (const int* it = lo; it != hi; ++it) {
                        int kk = static_cast<int>(it-mm_data);
                        out[*it - site_start] += Wi[miv[kk]] * my_signs[kk];
                    }
                }
            } else {
                for (int j = fi; j < li; j++) {
                    int a = tree_ivl_start[j], b = tree_ivl_end[j];
                    // Carry logic runs over all mismatches in [a,b) for
                    // correct state, but only writes to out for in-range sites.
                    // Skip entirely if interval doesn't overlap range.
                    int carry = (j==0)?0:static_cast<int>(tree_carry[static_cast<size_t>(i)*ni+j-1]);
                    const int* lo = std::lower_bound(mm_data, mm_data+n_mm, a);
                    const int* hi_mm = std::lower_bound(mm_data, mm_data+n_mm, b);
                    if (b <= site_start || a >= site_end) {
                        // Interval outside range: still update carry state
                        for (const int* it = lo; it != hi_mm; ++it) carry ^= 1;
                        continue;
                    }
                    double wi = Wi[j]; int prev = a;
                    for (const int* it = lo; it != hi_mm; ++it) {
                        int ms = *it;
                        if (carry && ms > prev) {
                            int rs = std::max(prev, site_start), re = std::min(ms, site_end);
                            for (int s = rs; s < re; s++) out[s-site_start] += wi;
                        }
                        carry ^= 1;
                        if (carry && ms >= site_start && ms < site_end)
                            out[ms-site_start] += wi;
                        prev = ms+1;
                    }
                    if (carry && b > prev) {
                        int rs = std::max(prev, site_start), re = std::min(b, site_end);
                        for (int s = rs; s < re; s++) out[s-site_start] += wi;
                    }
                }
            }
        }
    }

    // Sample 0
    { const double* W0 = &W[0];
      const int* ref = tree_ref_genome.data();
      for (int j = 0; j < ni; j++) {
          double w0j = W0[j];
          int a = std::max(tree_ivl_start[j], site_start);
          int b = std::min(tree_ivl_end[j], site_end);
          for (int s = a; s < b; s++) out[s-site_start] += w0j * ref[s];
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
