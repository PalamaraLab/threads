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
#include <cmath>
#include <cstring>
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
