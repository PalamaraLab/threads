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
#include <unordered_set>
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
    // target changes). Use the same chain-propagation as GenotypeIterator.
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
// of 1-run counts across all samples.

void ThreadingInstructions::prepare_rle_multiply() {
    if (rle_ready) return;

    const int n = num_samples;
    const int m = num_sites;

    rle_offset.resize(n + 1);
    std::vector<int> run_starts, run_ends;

    // Reference counting for genotype cache (shared with DAG prepare)
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

void ThreadingInstructions::compute_allele_counts() {
    if (allele_counts_ready) return;
    allele_counts.assign(num_sites, 0);
    GenotypeIterator gi(*this);
    int s = 0;
    while (gi.has_next_genotype()) {
        const std::vector<int>& g = gi.next_genotype();
        int ac = 0;
        for (int i = 0; i < num_samples; i++) ac += g[i];
        allele_counts[s++] = ac;
    }
    allele_counts_ready = true;
}

const std::vector<int>& ThreadingInstructions::get_allele_counts() {
    compute_allele_counts();
    return allele_counts;
}

// ---------------------------------------------------------------------------
// ARG traversal: visit_clades / visit_branches
// ---------------------------------------------------------------------------
//
// Reconstructs the implicit coalescent tree from threading instructions at each
// interval (contiguous region where no segment boundary changes the tree) and
// enumerates clades/branches via union-find.

namespace {

// Dendrogram node: O(n) memory representation of an implicit coalescent tree.
// n leaves (implicit, IDs 0..n-1) + up to n-1 internal nodes (IDs n..2n-2).
struct DNode {
    int left, right;    // child node indices
    double height;      // coalescence height
    uint64_t hash;      // commutative bottom-up hash for clade/branch identity
};

inline uint64_t dendro_leaf_hash(int id) {
    return static_cast<uint64_t>(id) * 0x9e3779b97f4a7c15ULL ^ 0x517cc1b727220a95ULL;
}

inline uint64_t dendro_combine(uint64_t h1, uint64_t h2, double height) {
    if (h1 > h2) std::swap(h1, h2);
    uint64_t hbits;
    std::memcpy(&hbits, &height, sizeof(double));
    uint64_t h = h1 * 0x9e3779b97f4a7c15ULL + h2;
    h ^= hbits + 0x517cc1b727220a95ULL + (h << 6) + (h >> 2);
    return h;
}

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

    // Emit leaf clades spanning the full region
    for (int i = 0; i < n; i++)
        callback({i}, 0.0, start, end);

    if (n <= 1) return;
    const int ne = n - 1;

    // Precompute leaf hashes (immutable, O(n))
    std::vector<uint64_t> lhash(n);
    for (int i = 0; i < n; i++)
        lhash[i] = dendro_leaf_hash(i);

    // Precompute edge events sorted by genomic position
    struct EdgeEvent { int bp; int sample; int target; double tmrca; };
    std::vector<EdgeEvent> events;
    events.reserve(ne * 20);
    for (int i = 1; i < n; i++) {
        const auto& instr = instructions[i];
        int seg0 = segment_at_bp(instr, start);
        events.push_back({start, i,
            seg0 < 0 ? 0 : instr.targets[seg0],
            seg0 < 0 ? 0.0 : instr.tmrcas[seg0]});
        for (int j = 0; j < static_cast<int>(instr.num_segments); j++) {
            int s = instr.starts[j];
            if (s <= start) continue;
            if (s >= end) break;
            events.push_back({s, i, instr.targets[j], instr.tmrcas[j]});
        }
    }
    std::sort(events.begin(), events.end(),
              [](const EdgeEvent& a, const EdgeEvent& b) { return a.bp < b.bp; });

    // Single dendrogram buffer + previous snapshot (avoids double-buffer swap)
    std::vector<DNode> dendro(ne), prev_dendro(ne);
    int dcount = 0, prev_dcount = 0;

    // UF arrays (reused per rebuild, O(n))
    std::vector<int> uf_p(n), uf_r(n), uf_dendro(n);

    // Current tree state
    std::vector<int> cur_tgt(n, 0);
    std::vector<double> cur_tmr(n, 0.0);

    // Compact-and-merge sorted edges: O(n + k log k) per boundary instead of O(n log n)
    struct SE { double h; int c, p; };
    std::vector<SE> sorted(ne), merge_buf(ne);
    std::vector<int> spos(n, -1);   // sample → position in sorted array
    std::vector<SE> new_entries;
    std::vector<int> changed;
    changed.reserve(32);

    // Pending clades: hash → start_bp. Pre-allocated.
    std::unordered_map<uint64_t, int> pending;
    pending.reserve(ne * 2);

    // Diff map: hash → old dendrogram node index. Pre-allocated, reused.
    std::unordered_map<uint64_t, int> old_map;
    old_map.reserve(ne * 2);

    // DFS buffers for on-demand descendant materialization (reused, O(n))
    std::vector<int> dfs_stack, desc_buf;
    dfs_stack.reserve(n);
    desc_buf.reserve(n);

    auto collect_desc = [&](const std::vector<DNode>& dn, int node_id) {
        desc_buf.clear();
        dfs_stack.clear();
        dfs_stack.push_back(node_id);
        while (!dfs_stack.empty()) {
            int id = dfs_stack.back();
            dfs_stack.pop_back();
            if (id < n) {
                desc_buf.push_back(id);
            } else {
                dfs_stack.push_back(dn[id - n].left);
                dfs_stack.push_back(dn[id - n].right);
            }
        }
        std::sort(desc_buf.begin(), desc_buf.end());
    };

    int ei = 0;
    const int E = static_cast<int>(events.size());
    bool need_init = true;
    bool have_prev = false;

    while (ei < E) {
        int bp = events[ei].bp;

        // Consume events, tracking changed samples
        changed.clear();
        bool tree_changed = false;
        while (ei < E && events[ei].bp == bp) {
            const auto& ev = events[ei++];
            if (need_init || ev.target != cur_tgt[ev.sample] || ev.tmrca != cur_tmr[ev.sample]) {
                cur_tgt[ev.sample] = ev.target;
                cur_tmr[ev.sample] = ev.tmrca;
                changed.push_back(ev.sample);
                tree_changed = true;
            }
        }
        if (!tree_changed) continue;

        // Compact-and-merge sorted edges
        if (need_init) {
            for (int i = 1; i < n; i++)
                sorted[i - 1] = {cur_tmr[i], i, cur_tgt[i]};
            std::sort(sorted.begin(), sorted.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            for (int i = 0; i < ne; i++) spos[sorted[i].c] = i;
            need_init = false;
        } else {
            std::sort(changed.begin(), changed.end(),
                      [&](int a, int b) { return spos[a] < spos[b]; });
            int blen = 0, ci = 0;
            for (int i = 0; i < ne; i++) {
                if (ci < static_cast<int>(changed.size()) && i == spos[changed[ci]])
                    ci++;
                else
                    merge_buf[blen++] = sorted[i];
            }
            new_entries.clear();
            for (int s : changed)
                new_entries.push_back({cur_tmr[s], s, cur_tgt[s]});
            std::sort(new_entries.begin(), new_entries.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            int ri = 0, ni = 0, out = 0;
            while (ri < blen && ni < static_cast<int>(new_entries.size())) {
                if (merge_buf[ri].h <= new_entries[ni].h)
                    sorted[out] = merge_buf[ri++];
                else
                    sorted[out] = new_entries[ni++];
                spos[sorted[out].c] = out;
                out++;
            }
            while (ri < blen) { sorted[out] = merge_buf[ri++]; spos[sorted[out].c] = out; out++; }
            while (ni < static_cast<int>(new_entries.size())) { sorted[out] = new_entries[ni++]; spos[sorted[out].c] = out; out++; }
        }

        // Build dendrogram via bottom-up UF merges
        for (int i = 0; i < n; i++) {
            uf_p[i] = i;
            uf_r[i] = 0;
            uf_dendro[i] = i;
        }
        dcount = 0;

        for (int i = 0; i < ne; i++) {
            int a = sorted[i].c, b = sorted[i].p;
            int ra = a;
            while (uf_p[ra] != ra) { uf_p[ra] = uf_p[uf_p[ra]]; ra = uf_p[ra]; }
            int rb = b;
            while (uf_p[rb] != rb) { uf_p[rb] = uf_p[uf_p[rb]]; rb = uf_p[rb]; }
            if (ra == rb) continue;

            int na = uf_dendro[ra], nb = uf_dendro[rb];
            uint64_t ha = (na < n) ? lhash[na] : dendro[na - n].hash;
            uint64_t hb = (nb < n) ? lhash[nb] : dendro[nb - n].hash;

            int inode = dcount++;
            dendro[inode] = {na, nb, sorted[i].h,
                             dendro_combine(ha, hb, sorted[i].h)};

            if (uf_r[ra] < uf_r[rb]) std::swap(ra, rb);
            uf_p[rb] = ra;
            if (uf_r[ra] == uf_r[rb]) uf_r[ra]++;
            uf_dendro[ra] = n + inode;
        }

        // Single-map diff: build map from old hashes, probe with new hashes
        if (have_prev) {
            old_map.clear();
            for (int i = 0; i < prev_dcount; i++)
                old_map[prev_dendro[i].hash] = i;

            for (int i = 0; i < dcount; i++) {
                uint64_t h = dendro[i].hash;
                auto it = old_map.find(h);
                if (it != old_map.end())
                    old_map.erase(it);   // persisted clade
                else
                    pending[h] = bp;     // new clade
            }

            // Flush disappeared clades (remaining in old_map)
            for (const auto& [h, node_idx] : old_map) {
                auto pit = pending.find(h);
                if (pit != pending.end()) {
                    collect_desc(prev_dendro, n + node_idx);
                    callback(desc_buf, prev_dendro[node_idx].height, pit->second, bp);
                    pending.erase(pit);
                }
            }
        } else {
            // First tree: all clades are new
            for (int i = 0; i < dcount; i++)
                pending[dendro[i].hash] = bp;
        }

        // Save current dendrogram for next iteration's diff
        std::copy(dendro.begin(), dendro.begin() + dcount, prev_dendro.begin());
        prev_dcount = dcount;
        have_prev = true;
    }

    // Flush all remaining pending clades
    for (int i = 0; i < prev_dcount; i++) {
        auto it = pending.find(prev_dendro[i].hash);
        if (it != pending.end()) {
            collect_desc(prev_dendro, n + i);
            callback(desc_buf, prev_dendro[i].height, it->second, end);
            pending.erase(it);
        }
    }
}


void ThreadingInstructions::visit_branches(
    std::function<void(const std::vector<int>&, double, double, int, int)> callback) const
{
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return;
    if (n <= 1) return;
    const int ne = n - 1;

    // Precompute leaf hashes
    std::vector<uint64_t> lhash(n);
    for (int i = 0; i < n; i++)
        lhash[i] = dendro_leaf_hash(i);

    // Precompute edge events sorted by genomic position
    struct EdgeEvent { int bp; int sample; int target; double tmrca; };
    std::vector<EdgeEvent> events;
    events.reserve(ne * 20);
    for (int i = 1; i < n; i++) {
        const auto& instr = instructions[i];
        int seg0 = segment_at_bp(instr, start);
        events.push_back({start, i,
            seg0 < 0 ? 0 : instr.targets[seg0],
            seg0 < 0 ? 0.0 : instr.tmrcas[seg0]});
        for (int j = 0; j < static_cast<int>(instr.num_segments); j++) {
            int s = instr.starts[j];
            if (s <= start) continue;
            if (s >= end) break;
            events.push_back({s, i, instr.targets[j], instr.tmrcas[j]});
        }
    }
    std::sort(events.begin(), events.end(),
              [](const EdgeEvent& a, const EdgeEvent& b) { return a.bp < b.bp; });

    // Single dendrogram buffer + previous snapshot
    std::vector<DNode> dendro(ne), prev_dendro(ne);
    int dcount = 0, prev_dcount = 0;

    // UF arrays
    std::vector<int> uf_p(n), uf_r(n), uf_dendro(n);

    // Current tree state
    std::vector<int> cur_tgt(n, 0);
    std::vector<double> cur_tmr(n, 0.0);

    // Compact-and-merge sorted edges
    struct SE { double h; int c, p; };
    std::vector<SE> sorted(ne), merge_buf(ne);
    std::vector<int> spos(n, -1);
    std::vector<SE> new_entries;
    std::vector<int> changed;
    changed.reserve(32);

    // Branch hash
    auto branch_hash_fn = [](uint64_t child_hash, double child_h, double parent_h) -> uint64_t {
        uint64_t ch_bits, ph_bits;
        std::memcpy(&ch_bits, &child_h, sizeof(double));
        std::memcpy(&ph_bits, &parent_h, sizeof(double));
        uint64_t h = child_hash ^ (ch_bits * 0x9e3779b97f4a7c15ULL);
        h ^= ph_bits + 0x517cc1b727220a95ULL + (h << 6) + (h >> 2);
        return h;
    };

    // Branch info for enumeration from dendrogram
    struct BranchInfo {
        uint64_t hash;
        int child_node;
        double child_h, parent_h;
    };
    std::vector<BranchInfo> branches, prev_branches;
    branches.reserve(2 * ne);
    prev_branches.reserve(2 * ne);

    auto enumerate_branches = [&](const std::vector<DNode>& dn, int count,
                                   std::vector<BranchInfo>& out) {
        out.clear();
        for (int i = 0; i < count; i++) {
            const DNode& nd = dn[i];
            double lh = (nd.left < n) ? 0.0 : dn[nd.left - n].height;
            if (lh < nd.height) {
                uint64_t ch = (nd.left < n) ? lhash[nd.left] : dn[nd.left - n].hash;
                out.push_back({branch_hash_fn(ch, lh, nd.height), nd.left, lh, nd.height});
            }
            double rh = (nd.right < n) ? 0.0 : dn[nd.right - n].height;
            if (rh < nd.height) {
                uint64_t ch = (nd.right < n) ? lhash[nd.right] : dn[nd.right - n].hash;
                out.push_back({branch_hash_fn(ch, rh, nd.height), nd.right, rh, nd.height});
            }
        }
    };

    // Pending branches: hash → start_bp. Pre-allocated.
    std::unordered_map<uint64_t, int> pending;
    pending.reserve(2 * ne * 2);

    // Diff map: pre-allocated, reused
    std::unordered_map<uint64_t, int> old_map;  // hash → index in prev_branches
    old_map.reserve(2 * ne * 2);

    // DFS buffers
    std::vector<int> dfs_stack, desc_buf;
    dfs_stack.reserve(n);
    desc_buf.reserve(n);

    auto collect_desc = [&](const std::vector<DNode>& dn, int node_id) {
        desc_buf.clear();
        dfs_stack.clear();
        dfs_stack.push_back(node_id);
        while (!dfs_stack.empty()) {
            int id = dfs_stack.back();
            dfs_stack.pop_back();
            if (id < n) {
                desc_buf.push_back(id);
            } else {
                dfs_stack.push_back(dn[id - n].left);
                dfs_stack.push_back(dn[id - n].right);
            }
        }
        std::sort(desc_buf.begin(), desc_buf.end());
    };

    int ei = 0;
    const int E = static_cast<int>(events.size());
    bool need_init = true;
    bool have_prev = false;

    while (ei < E) {
        int bp = events[ei].bp;

        changed.clear();
        bool tree_changed = false;
        while (ei < E && events[ei].bp == bp) {
            const auto& ev = events[ei++];
            if (need_init || ev.target != cur_tgt[ev.sample] || ev.tmrca != cur_tmr[ev.sample]) {
                cur_tgt[ev.sample] = ev.target;
                cur_tmr[ev.sample] = ev.tmrca;
                changed.push_back(ev.sample);
                tree_changed = true;
            }
        }
        if (!tree_changed) continue;

        // Compact-and-merge sorted edges
        if (need_init) {
            for (int i = 1; i < n; i++)
                sorted[i - 1] = {cur_tmr[i], i, cur_tgt[i]};
            std::sort(sorted.begin(), sorted.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            for (int i = 0; i < ne; i++) spos[sorted[i].c] = i;
            need_init = false;
        } else {
            std::sort(changed.begin(), changed.end(),
                      [&](int a, int b) { return spos[a] < spos[b]; });
            int blen = 0, ci = 0;
            for (int i = 0; i < ne; i++) {
                if (ci < static_cast<int>(changed.size()) && i == spos[changed[ci]])
                    ci++;
                else
                    merge_buf[blen++] = sorted[i];
            }
            new_entries.clear();
            for (int s : changed)
                new_entries.push_back({cur_tmr[s], s, cur_tgt[s]});
            std::sort(new_entries.begin(), new_entries.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            int ri = 0, ni = 0, out = 0;
            while (ri < blen && ni < static_cast<int>(new_entries.size())) {
                if (merge_buf[ri].h <= new_entries[ni].h)
                    sorted[out] = merge_buf[ri++];
                else
                    sorted[out] = new_entries[ni++];
                spos[sorted[out].c] = out;
                out++;
            }
            while (ri < blen) { sorted[out] = merge_buf[ri++]; spos[sorted[out].c] = out; out++; }
            while (ni < static_cast<int>(new_entries.size())) { sorted[out] = new_entries[ni++]; spos[sorted[out].c] = out; out++; }
        }

        // Build dendrogram via UF
        for (int i = 0; i < n; i++) {
            uf_p[i] = i;
            uf_r[i] = 0;
            uf_dendro[i] = i;
        }
        dcount = 0;

        for (int i = 0; i < ne; i++) {
            int a = sorted[i].c, b = sorted[i].p;
            int ra = a;
            while (uf_p[ra] != ra) { uf_p[ra] = uf_p[uf_p[ra]]; ra = uf_p[ra]; }
            int rb = b;
            while (uf_p[rb] != rb) { uf_p[rb] = uf_p[uf_p[rb]]; rb = uf_p[rb]; }
            if (ra == rb) continue;

            int na = uf_dendro[ra], nb = uf_dendro[rb];
            uint64_t ha = (na < n) ? lhash[na] : dendro[na - n].hash;
            uint64_t hb = (nb < n) ? lhash[nb] : dendro[nb - n].hash;

            int inode = dcount++;
            dendro[inode] = {na, nb, sorted[i].h,
                             dendro_combine(ha, hb, sorted[i].h)};

            if (uf_r[ra] < uf_r[rb]) std::swap(ra, rb);
            uf_p[rb] = ra;
            if (uf_r[ra] == uf_r[rb]) uf_r[ra]++;
            uf_dendro[ra] = n + inode;
        }

        // Single-map diff for branches
        if (have_prev) {
            old_map.clear();
            for (int i = 0; i < static_cast<int>(prev_branches.size()); i++)
                old_map[prev_branches[i].hash] = i;

            enumerate_branches(dendro, dcount, branches);
            for (const auto& br : branches) {
                auto it = old_map.find(br.hash);
                if (it != old_map.end())
                    old_map.erase(it);   // persisted branch
                else
                    pending[br.hash] = bp;  // new branch
            }

            // Flush disappeared branches
            for (const auto& [h, idx] : old_map) {
                auto pit = pending.find(h);
                if (pit != pending.end()) {
                    const auto& br = prev_branches[idx];
                    collect_desc(prev_dendro, br.child_node);
                    callback(desc_buf, br.child_h, br.parent_h, pit->second, bp);
                    pending.erase(pit);
                }
            }
        } else {
            enumerate_branches(dendro, dcount, branches);
            for (const auto& br : branches)
                pending[br.hash] = bp;
        }

        // Save current state for next iteration
        std::copy(dendro.begin(), dendro.begin() + dcount, prev_dendro.begin());
        prev_dcount = dcount;
        prev_branches.assign(branches.begin(), branches.end());
        have_prev = true;
    }

    // Flush remaining branches
    if (have_prev) {
        for (const auto& br : prev_branches) {
            auto it = pending.find(br.hash);
            if (it != pending.end()) {
                collect_desc(prev_dendro, br.child_node);
                callback(desc_buf, br.child_h, br.parent_h, it->second, end);
                pending.erase(it);
            }
        }
    }
}




double ThreadingInstructions::total_volume() const {
    const int n = num_samples;
    if (n <= 1 || num_sites == 0) return 0.0;

    // Key insight: total branch length of a coalescent tree built from n-1 merges
    // at sorted heights h_1 <= ... <= h_{n-1} equals:
    //
    //   L = sum_{j=0}^{n-2} (h_{j+1} - h_j) * (n - j) = sum(h) + max(h)
    //
    // (telescoping sum). Since threading edges form a tree (target[i] < i), every
    // edge creates a unique non-redundant merge (trees have no cycles, so adding
    // any tree edge to a forest of tree edges always connects two components).
    // Therefore the merge heights are exactly the n-1 tmrca values, and:
    //
    //   vol(position) = sum(tmrcas) + max(tmrcas)
    //
    // No UF or sorting needed — just a running sum and max, updated in O(k) when
    // k edges change at a boundary. Total work: O(E) for E total edge events.

    // Precompute edge-change events sorted by position.
    struct EdgeEvent { int bp; int sample; double tmrca; };
    std::vector<EdgeEvent> events;
    events.reserve((n - 1) * 20);

    for (int i = 1; i < n; i++) {
        const auto& instr = instructions[i];
        int seg0 = segment_at_bp(instr, start);
        events.push_back({start, i, seg0 < 0 ? 0.0 : instr.tmrcas[seg0]});
        for (int j = 0; j < static_cast<int>(instr.num_segments); j++) {
            int s = instr.starts[j];
            if (s <= start) continue;
            if (s >= end) break;
            events.push_back({s, i, instr.tmrcas[j]});
        }
    }
    std::sort(events.begin(), events.end(),
              [](const EdgeEvent& a, const EdgeEvent& b) { return a.bp < b.bp; });

    // Sweep: maintain sum and max of tmrcas incrementally.
    std::vector<double> cur_h(n, 0.0);
    double sum_h = 0.0;
    double max_h = 0.0;

    double vol = 0.0;
    int ei = 0;
    const int E = static_cast<int>(events.size());

    while (ei < E) {
        int bp = events[ei].bp;
        bool max_possibly_removed = false;
        double new_max_candidate = 0.0;

        // Apply all events at this position
        while (ei < E && events[ei].bp == bp) {
            int s = events[ei].sample;
            double new_h = events[ei].tmrca;
            double old_h = cur_h[s];
            ei++;
            if (new_h != old_h) {
                sum_h += (new_h - old_h);
                if (old_h >= max_h) max_possibly_removed = true;
                cur_h[s] = new_h;
                if (new_h > new_max_candidate) new_max_candidate = new_h;
            }
        }

        // Fix max: O(n) rescan only when current max was removed; O(1) otherwise
        if (max_possibly_removed) {
            max_h = *std::max_element(cur_h.begin() + 1, cur_h.end());
        } else if (new_max_candidate > max_h) {
            max_h = new_max_candidate;
        }

        int bp_next = (ei < E) ? events[ei].bp : end;
        double span = static_cast<double>(bp_next - bp);
        vol += (sum_h + max_h) * span;
    }
    return vol;
}


std::vector<double> ThreadingInstructions::allele_frequency_spectrum() const {
    const int n = num_samples;
    if (n == 0 || num_sites == 0) return std::vector<double>(n + 1, 0.0);
    if (n == 1) return std::vector<double>(2, 0.0);

    const int ne = n - 1;

    // Precompute edge-change events (same strategy as total_volume)
    struct EdgeEvent { int bp; int sample; int target; double tmrca; };
    std::vector<EdgeEvent> events;
    events.reserve(ne * 20);

    for (int i = 1; i < n; i++) {
        const auto& instr = instructions[i];
        int seg0 = segment_at_bp(instr, start);
        events.push_back({start, i,
            seg0 < 0 ? 0 : instr.targets[seg0],
            seg0 < 0 ? 0.0 : instr.tmrcas[seg0]});
        for (int j = 0; j < static_cast<int>(instr.num_segments); j++) {
            int s = instr.starts[j];
            if (s <= start) continue;
            if (s >= end) break;
            events.push_back({s, i, instr.targets[j], instr.tmrcas[j]});
        }
    }
    std::sort(events.begin(), events.end(),
              [](const EdgeEvent& a, const EdgeEvent& b) { return a.bp < b.bp; });

    std::vector<int> cur_tgt(n, 0);
    std::vector<double> cur_tmr(n, 0.0);

    struct SE { double h; int c, p; };
    std::vector<SE> sorted(ne), buf(ne);
    std::vector<int> spos(n, -1);

    // UF with sizes for frequency spectrum binning
    std::vector<int> uf_p(n), uf_r(n), uf_sz(n);
    std::vector<double> uf_h(n);
    auto uf_find = [&](int x) {
        while (uf_p[x] != x) { uf_p[x] = uf_p[uf_p[x]]; x = uf_p[x]; }
        return x;
    };

    std::vector<SE> new_entries;
    std::vector<int> changed;
    changed.reserve(32);

    std::vector<double> afs(n + 1, 0.0);
    std::vector<double> prev_afs(n + 1, 0.0);
    // Track which bins are non-zero for fast span multiplication
    std::vector<int> active_bins;
    active_bins.reserve(n);

    int ei = 0;
    const int E = static_cast<int>(events.size());
    bool need_init = true;

    while (ei < E) {
        int bp = events[ei].bp;

        changed.clear();
        bool tree_changed = false;
        while (ei < E && events[ei].bp == bp) {
            const auto& ev = events[ei++];
            if (need_init || ev.target != cur_tgt[ev.sample] || ev.tmrca != cur_tmr[ev.sample]) {
                cur_tgt[ev.sample] = ev.target;
                cur_tmr[ev.sample] = ev.tmrca;
                changed.push_back(ev.sample);
                tree_changed = true;
            }
        }

        int bp_next = (ei < E) ? events[ei].bp : end;
        double span = static_cast<double>(bp_next - bp);

        if (!tree_changed) {
            for (int k : active_bins) afs[k] += prev_afs[k] * span;
            continue;
        }

        if (need_init) {
            for (int i = 1; i < n; i++)
                sorted[i - 1] = {cur_tmr[i], i, cur_tgt[i]};
            std::sort(sorted.begin(), sorted.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            for (int i = 0; i < ne; i++) spos[sorted[i].c] = i;
            need_init = false;
        } else {
            std::sort(changed.begin(), changed.end(),
                      [&](int a, int b) { return spos[a] < spos[b]; });
            int blen = 0, ci = 0;
            for (int i = 0; i < ne; i++) {
                if (ci < static_cast<int>(changed.size()) && i == spos[changed[ci]]) {
                    ci++;
                } else {
                    buf[blen++] = sorted[i];
                }
            }
            new_entries.clear();
            for (int s : changed)
                new_entries.push_back({cur_tmr[s], s, cur_tgt[s]});
            std::sort(new_entries.begin(), new_entries.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            int ri = 0, ni = 0, out = 0;
            while (ri < blen && ni < static_cast<int>(new_entries.size())) {
                if (buf[ri].h <= new_entries[ni].h) {
                    sorted[out] = buf[ri++];
                } else {
                    sorted[out] = new_entries[ni++];
                }
                spos[sorted[out].c] = out;
                out++;
            }
            while (ri < blen) { sorted[out] = buf[ri++]; spos[sorted[out].c] = out; out++; }
            while (ni < static_cast<int>(new_entries.size())) { sorted[out] = new_entries[ni++]; spos[sorted[out].c] = out; out++; }
        }

        // UF rebuild with size tracking
        for (int i = 0; i < n; i++) { uf_p[i] = i; uf_r[i] = 0; uf_sz[i] = 1; uf_h[i] = 0.0; }

        // Clear previous AFS contributions (only non-zero bins)
        for (int k : active_bins) prev_afs[k] = 0.0;
        active_bins.clear();

        for (int i = 0; i < ne; i++) {
            int a = uf_find(sorted[i].c), b = uf_find(sorted[i].p);
            if (a != b) {
                int sa = uf_sz[a], sb = uf_sz[b];
                double ha = uf_h[a], hb = uf_h[b];
                if (ha < sorted[i].h) prev_afs[sa] += sorted[i].h - ha;
                if (hb < sorted[i].h) prev_afs[sb] += sorted[i].h - hb;
                if (uf_r[a] < uf_r[b]) std::swap(a, b);
                uf_p[b] = a;
                uf_sz[a] += uf_sz[b];
                if (uf_r[a] == uf_r[b]) uf_r[a]++;
                uf_h[a] = sorted[i].h;
            }
        }

        // Collect active bins and accumulate
        for (int k = 1; k < n; k++) {
            if (prev_afs[k] > 0.0) {
                active_bins.push_back(k);
                afs[k] += prev_afs[k] * span;
            }
        }
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
    const int ne = n - 1;

    // Event-based boundary processing (same pattern as AFS).
    // Key optimization vs visit_branches: no hash maps, no persistence tracking,
    // each boundary is processed independently. Only 1 Poisson sample per boundary
    // (thinning), and DFS only for the rare branches that actually get mutations.
    struct EdgeEvent { int bp; int sample; int target; double tmrca; };
    std::vector<EdgeEvent> events;
    events.reserve(ne * 20);
    for (int i = 1; i < n; i++) {
        const auto& instr = instructions[i];
        int seg0 = segment_at_bp(instr, start);
        events.push_back({start, i,
            seg0 < 0 ? 0 : instr.targets[seg0],
            seg0 < 0 ? 0.0 : instr.tmrcas[seg0]});
        for (int j = 0; j < static_cast<int>(instr.num_segments); j++) {
            int s = instr.starts[j];
            if (s <= start) continue;
            if (s >= end) break;
            events.push_back({s, i, instr.targets[j], instr.tmrcas[j]});
        }
    }
    std::sort(events.begin(), events.end(),
              [](const EdgeEvent& a, const EdgeEvent& b) { return a.bp < b.bp; });

    std::vector<int> cur_tgt(n, 0);
    std::vector<double> cur_tmr(n, 0.0);

    // Compact-and-merge sorted edges
    struct SE { double h; int c, p; };
    std::vector<SE> sorted(ne), merge_buf(ne);
    std::vector<int> spos(n, -1);
    std::vector<SE> new_entries;
    std::vector<int> changed;
    changed.reserve(32);

    // UF arrays with sizes and heights (no hashing needed)
    std::vector<int> uf_p(n), uf_r(n), uf_dendro(n);
    std::vector<double> uf_h(n);
    auto uf_find = [&](int x) {
        while (uf_p[x] != x) { uf_p[x] = uf_p[uf_p[x]]; x = uf_p[x]; }
        return x;
    };

    // Dendrogram for DFS (built alongside UF, no hash computation)
    std::vector<DNode> dendro(ne);

    // Branch records for weighted Poisson thinning
    struct BranchRecord { int node; double vol; };
    std::vector<BranchRecord> br_rec;
    br_rec.reserve(2 * ne);
    std::vector<double> cum_vol;
    cum_vol.reserve(2 * ne);

    // DFS buffers
    std::vector<int> dfs_stack, desc_buf;
    dfs_stack.reserve(n);
    desc_buf.reserve(n);

    int ei = 0;
    const int E = static_cast<int>(events.size());
    bool need_init = true;

    while (ei < E) {
        int bp = events[ei].bp;

        changed.clear();
        bool tree_changed = false;
        while (ei < E && events[ei].bp == bp) {
            const auto& ev = events[ei++];
            if (need_init || ev.target != cur_tgt[ev.sample] || ev.tmrca != cur_tmr[ev.sample]) {
                cur_tgt[ev.sample] = ev.target;
                cur_tmr[ev.sample] = ev.tmrca;
                changed.push_back(ev.sample);
                tree_changed = true;
            }
        }

        int bp_next = (ei < E) ? events[ei].bp : end;
        double span = static_cast<double>(bp_next - bp);
        if (!tree_changed || span <= 0) continue;

        // Compact-and-merge sorted edges
        if (need_init) {
            for (int i = 1; i < n; i++)
                sorted[i - 1] = {cur_tmr[i], i, cur_tgt[i]};
            std::sort(sorted.begin(), sorted.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            for (int i = 0; i < ne; i++) spos[sorted[i].c] = i;
            need_init = false;
        } else {
            std::sort(changed.begin(), changed.end(),
                      [&](int a, int b) { return spos[a] < spos[b]; });
            int blen = 0, ci = 0;
            for (int i = 0; i < ne; i++) {
                if (ci < static_cast<int>(changed.size()) && i == spos[changed[ci]])
                    ci++;
                else
                    merge_buf[blen++] = sorted[i];
            }
            new_entries.clear();
            for (int s : changed)
                new_entries.push_back({cur_tmr[s], s, cur_tgt[s]});
            std::sort(new_entries.begin(), new_entries.end(),
                      [](const SE& a, const SE& b) { return a.h < b.h; });
            int ri = 0, ni = 0, out = 0;
            while (ri < blen && ni < static_cast<int>(new_entries.size())) {
                if (merge_buf[ri].h <= new_entries[ni].h)
                    sorted[out] = merge_buf[ri++];
                else
                    sorted[out] = new_entries[ni++];
                spos[sorted[out].c] = out;
                out++;
            }
            while (ri < blen) { sorted[out] = merge_buf[ri++]; spos[sorted[out].c] = out; out++; }
            while (ni < static_cast<int>(new_entries.size())) { sorted[out] = new_entries[ni++]; spos[sorted[out].c] = out; out++; }
        }

        // UF with dendrogram build + branch volume collection
        for (int i = 0; i < n; i++) {
            uf_p[i] = i; uf_r[i] = 0; uf_h[i] = 0.0; uf_dendro[i] = i;
        }
        int dcount = 0;
        br_rec.clear();

        for (int i = 0; i < ne; i++) {
            int a = uf_find(sorted[i].c), b = uf_find(sorted[i].p);
            if (a == b) continue;

            double h = sorted[i].h;
            double ha = uf_h[a], hb = uf_h[b];
            int na = uf_dendro[a], nb = uf_dendro[b];

            // Build dendrogram node (no hash needed)
            dendro[dcount] = {na, nb, h, 0};

            // Collect branch volumes for Poisson thinning
            double vol_a = (h - ha) * span;
            if (vol_a > 0) br_rec.push_back({na, vol_a});
            double vol_b = (h - hb) * span;
            if (vol_b > 0) br_rec.push_back({nb, vol_b});

            // UF merge
            if (uf_r[a] < uf_r[b]) std::swap(a, b);
            uf_p[b] = a;
            if (uf_r[a] == uf_r[b]) uf_r[a]++;
            uf_h[a] = h;
            uf_dendro[a] = n + dcount;
            dcount++;
        }

        if (br_rec.empty()) continue;

        // Poisson thinning: one sample for total volume, then distribute
        double total_vol = 0;
        cum_vol.resize(br_rec.size());
        for (size_t i = 0; i < br_rec.size(); i++) {
            total_vol += br_rec[i].vol;
            cum_vol[i] = total_vol;
        }

        double expected = mutation_rate * total_vol;
        if (expected <= 0) continue;

        std::poisson_distribution<int> pois(expected);
        int n_muts = pois(rng);

        for (int m = 0; m < n_muts; m++) {
            // Choose branch proportional to volume
            double u = std::uniform_real_distribution<>(0, total_vol)(rng);
            int idx = static_cast<int>(
                std::lower_bound(cum_vol.begin(), cum_vol.end(), u) - cum_vol.begin());
            if (idx >= static_cast<int>(br_rec.size()))
                idx = static_cast<int>(br_rec.size()) - 1;

            // Choose position within span
            int pos = bp + static_cast<int>(
                std::uniform_real_distribution<double>(0, span)(rng));

            // DFS for descendants (only for branches that get mutations)
            desc_buf.clear();
            dfs_stack.clear();
            dfs_stack.push_back(br_rec[idx].node);
            while (!dfs_stack.empty()) {
                int id = dfs_stack.back();
                dfs_stack.pop_back();
                if (id < n) {
                    desc_buf.push_back(id);
                } else {
                    dfs_stack.push_back(dendro[id - n].left);
                    dfs_stack.push_back(dendro[id - n].right);
                }
            }

            // Add mutation to result
            int offset = result.n_mutations * n;
            result.genotypes.resize(offset + n, 0);
            for (int d : desc_buf)
                result.genotypes[offset + d] = 1;
            result.positions.push_back(pos);
            result.n_mutations++;
        }
    }

    return result;
}


void ThreadingInstructions::visit_mutations(
    std::function<void(int, const std::vector<int>&)> callback) const
{
    if (num_samples == 0 || num_sites == 0) return;

    GenotypeIterator gi(*this);
    std::vector<int> carriers;

    for (int s = 0; s < num_sites; s++) {
        const std::vector<int>& g = gi.next_genotype();
        carriers.clear();
        for (int i = 0; i < num_samples; i++) {
            if (g[i] == 1) carriers.push_back(i);
        }
        if (!carriers.empty()) {
            callback(positions[s], carriers);
        }
    }
}


double ThreadingInstructions::mrca(int id1, int id2, int position_bp) const {
    if (id1 < 0 || id1 >= num_samples || id2 < 0 || id2 >= num_samples) {
        throw std::runtime_error("Sample ID out of range");
    }
    if (id1 == id2) return 0.0;

    // Trace path from id1 to root, recording max coalescence height at each node
    std::unordered_map<int, double> path1;
    path1.reserve(64);
    double h1 = 0.0;
    int cur = id1;
    path1[cur] = 0.0;
    while (cur != 0) {
        int seg = segment_at_bp(instructions[cur], position_bp);
        int tgt = (seg < 0) ? 0 : instructions[cur].targets[seg];
        double edge_h = (seg < 0) ? 0.0 : instructions[cur].tmrcas[seg];
        h1 = std::max(h1, edge_h);
        cur = tgt;
        path1[cur] = h1;
    }

    // Walk from id2 upward until hitting a node on path1
    double h2 = 0.0;
    cur = id2;
    int steps = 0;
    while (path1.find(cur) == path1.end()) {
        int seg = segment_at_bp(instructions[cur], position_bp);
        int tgt = (seg < 0) ? 0 : instructions[cur].targets[seg];
        double edge_h = (seg < 0) ? 0.0 : instructions[cur].tmrcas[seg];
        h2 = std::max(h2, edge_h);
        cur = tgt;
        if (++steps > num_samples) {
            // Safety: avoid infinite loop on malformed instructions
            return std::max(h1, h2);
        }
    }

    return std::max(path1[cur], h2);
}


std::vector<double> ThreadingInstructions::mrca_vector(int id1, int id2) const {
    std::vector<double> result(num_sites);
    for (int s = 0; s < num_sites; s++) {
        result[s] = mrca(id1, id2, positions[s]);
    }
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
