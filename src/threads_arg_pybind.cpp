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

#include "ImputationMatcher.hpp"
#include "ThreadsLowMem.hpp"
#include "DataConsistency.hpp"
#include "AlleleAges.hpp"
#include "GenotypeIterator.hpp"
#include "VCFWriter.hpp"
#include "ForwardBackward.hpp"
#include "ThreadsIO.hpp"
#include "pybind_utils.hpp"

#include <pybind11/numpy.h>
#include <cstring>
#include <vector>

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Normalize/diploid helpers for DAG and RLE multiply wrappers.
// All pre/post-processing is O(m) or O(n), negligible vs the DAG traversal.
// ---------------------------------------------------------------------------

struct NormStats {
    std::vector<double> inv_std;  // 1/std per site
    std::vector<double> mean;     // mean per site (haploid p or diploid 2p)
};

// Compute per-site normalization stats.
// Haploid: mean = p, std = sqrt(p*(1-p))
// Diploid: mean = ac/n_dip, std = empirical dosage std (matches dense baseline)
static NormStats get_norm_stats(ThreadingInstructions& self, bool diploid) {
    self.compute_allele_counts();
    const auto& ac = self.get_allele_counts();
    const int n = self.num_samples;
    const int m = self.num_sites;
    NormStats ns;
    ns.inv_std.resize(m);
    ns.mean.resize(m);
    if (diploid) {
        const int n_dip = n / 2;
        // Need het counts to compute empirical dosage variance:
        // dosage_sumsq = 2*ac - n_het, var = sumsq/n_dip - mean^2
        auto het = self.het_per_site();
        for (int s = 0; s < m; s++) {
            double mu = static_cast<double>(ac[s]) / n_dip;
            double sumsq = static_cast<double>(2 * ac[s] - het[s]);
            double var = sumsq / n_dip - mu * mu;
            ns.mean[s] = mu;
            ns.inv_std[s] = (var > 0) ? 1.0 / std::sqrt(var) : 0.0;
        }
    } else {
        for (int s = 0; s < m; s++) {
            double p = static_cast<double>(ac[s]) / n;
            double var = p * (1.0 - p);
            ns.mean[s] = p;
            ns.inv_std[s] = (var > 0) ? 1.0 / std::sqrt(var) : 0.0;
        }
    }
    return ns;
}

// Pre-process x (sites) for normalized right multiply: x_scaled[s] = x[s] * inv_std[s]
// Also computes offset = sum_s(mean[s] * x_scaled[s]) to subtract from result.
// Works on flat arrays with k columns: x[s*k + c].
static void norm_prescale_right(std::vector<double>& x, const NormStats& ns,
                                std::vector<double>& offsets, int m, int k) {
    offsets.assign(k, 0.0);
    for (int s = 0; s < m; s++) {
        for (int c = 0; c < k; c++) {
            x[s * k + c] *= ns.inv_std[s];
            offsets[c] += ns.mean[s] * x[s * k + c];
        }
    }
}

// Post-process result of right multiply: reduce to diploid first, then subtract offset.
// result is flat (n, k). Output is (out_n, k) where out_n = n or n/2.
static void norm_postprocess_right(std::vector<double>& result, const std::vector<double>& offsets,
                                   bool diploid, int n, int k) {
    // Reduce to diploid first (offset must be subtracted once per diploid, not per haploid)
    int out_n = n;
    if (diploid) {
        const int n_dip = n / 2;
        for (int j = 0; j < n_dip; j++)
            for (int c = 0; c < k; c++)
                result[j * k + c] = result[(2*j) * k + c] + result[(2*j+1) * k + c];
        result.resize(static_cast<size_t>(n_dip) * k);
        out_n = n_dip;
    }
    // Subtract mean offset
    for (int i = 0; i < out_n; i++)
        for (int c = 0; c < k; c++)
            result[i * k + c] -= offsets[c];
}

// Reduce haploid result to diploid (no normalize).
static void diploid_reduce_right(std::vector<double>& result, int n, int k) {
    const int n_dip = n / 2;
    for (int j = 0; j < n_dip; j++)
        for (int c = 0; c < k; c++)
            result[j * k + c] = result[(2*j) * k + c] + result[(2*j+1) * k + c];
    result.resize(static_cast<size_t>(n_dip) * k);
}

// Pre-process y (samples) for left multiply with diploid: expand n_dip → n.
// y[j*k + c] → y_hap[(2j)*k + c] = y_hap[(2j+1)*k + c] = y[j*k + c]
static std::vector<double> diploid_expand_left(const std::vector<double>& y, int n, int k) {
    const int n_dip = n / 2;
    std::vector<double> y_hap(static_cast<size_t>(n) * k);
    for (int j = 0; j < n_dip; j++)
        for (int c = 0; c < k; c++) {
            y_hap[(2*j) * k + c] = y[j * k + c];
            y_hap[(2*j+1) * k + c] = y[j * k + c];
        }
    return y_hap;
}

// Post-process result of left multiply for normalize:
// result[s] = (result[s] - mean[s] * sum_y) / std[s]
static void norm_postprocess_left(std::vector<double>& result, const NormStats& ns,
                                  const std::vector<double>& col_sums, int m, int k) {
    for (int s = 0; s < m; s++)
        for (int c = 0; c < k; c++)
            result[s * k + c] = (result[s * k + c] - ns.mean[s] * col_sums[c]) * ns.inv_std[s];
}

// Compute column sums of input y (flat, rows × k).
static std::vector<double> col_sums_of(const double* y, int rows, int k) {
    std::vector<double> sums(k, 0.0);
    for (int i = 0; i < rows; i++)
        for (int c = 0; c < k; c++)
            sums[c] += y[i * k + c];
    return sums;
}

PYBIND11_MODULE(threads_arg_python_bindings, m) {
  py::class_<ImputationSegment>(m, "ImputationSegment")
      .def_readonly("seg_start", &ImputationSegment::seg_start)
      .def_readonly("ids", &ImputationSegment::ids)
      .def_readonly("weights", &ImputationSegment::weights)
      .def_readonly("ages", &ImputationSegment::ages);

  py::class_<ThreadsLowMem>(m, "ThreadsLowMem")
      .def(py::init<std::vector<int>, std::vector<double>, std::vector<double>, std::vector<double>,
                    std::vector<double>, double, bool>(),
           "Initialize", py::arg("target_ids"), py::arg("physical_positions"),
           py::arg("genetic_positions"), py::arg("ne_sizes"), py::arg("ne_times"),
           py::arg("mutation_rate"), py::arg("sparse"))
      .def_readonly("physical_positions", &ThreadsLowMem::physical_positions)
      .def_readonly("genetic_positions", &ThreadsLowMem::genetic_positions)
      .def_readonly("num_samples", &ThreadsLowMem::num_samples)
      .def_readonly("num_sites", &ThreadsLowMem::num_sites)
      .def_readonly("paths", &ThreadsLowMem::paths)
      .def_readonly("mean_bp_size", &ThreadsLowMem::mean_bp_size)
      .def_readonly("mutation_rate", &ThreadsLowMem::mutation_rate)
      .def_readonly("expected_branch_lengths", &ThreadsLowMem::expected_branch_lengths)
      .def("initialize_viterbi", &ThreadsLowMem::initialize_viterbi)
      .def("process_site_viterbi", &ThreadsLowMem::process_site_viterbi)
      .def("process_all_sites_viterbi", &ThreadsLowMem::process_all_sites_viterbi)
      .def("process_all_sites_viterbi_numpy", [](ThreadsLowMem& self, py::array_t<int32_t, py::array::c_style | py::array::forcecast> arr) {
        auto buf = arr.request();
        if (buf.ndim != 2) throw std::runtime_error("Expected 2D array (n_sites × n_haps)");
        int n_sites = static_cast<int>(buf.shape[0]);
        int n_haps = static_cast<int>(buf.shape[1]);
        self.process_all_sites_viterbi_flat(static_cast<const int32_t*>(buf.ptr), n_sites, n_haps);
      })
      .def("process_site_hets", &ThreadsLowMem::process_site_hets)
      .def("process_all_sites_hets", &ThreadsLowMem::process_all_sites_hets)
      .def("process_all_sites_hets_numpy", [](ThreadsLowMem& self, py::array_t<int32_t, py::array::c_style | py::array::forcecast> arr) {
        auto buf = arr.request();
        if (buf.ndim != 2) throw std::runtime_error("Expected 2D array (n_sites × n_haps)");
        int n_sites = static_cast<int>(buf.shape[0]);
        int n_haps = static_cast<int>(buf.shape[1]);
        self.process_all_sites_hets_flat(static_cast<const int32_t*>(buf.ptr), n_sites, n_haps);
      })
      .def("count_branches", &ThreadsLowMem::count_branches)
      .def("prune", &ThreadsLowMem::prune)
      .def("traceback", &ThreadsLowMem::traceback)
      .def("date_segments", &ThreadsLowMem::date_segments)
      .def("serialize_paths", &ThreadsLowMem::serialize_paths);

  py::class_<ViterbiPath>(m, "ViterbiPath")
      .def(py::init<int>(), "Initialize", py::arg("target_id"))
      .def(py::init<int, std::vector<int>, std::vector<int>, std::vector<double>,
                    std::vector<int>>(),
           "Initialize", py::arg("target_id"), py::arg("sample_ids"), py::arg("segment_starts"),
           py::arg("heights"), py::arg("het_sites"))
      .def_readonly("segment_starts", &ViterbiPath::segment_starts)
      .def_readonly("het_sites", &ViterbiPath::het_sites)
      .def_readonly("heights", &ViterbiPath::heights)
      .def_readonly("score", &ViterbiPath::score)
      .def_readonly("sample_ids", &ViterbiPath::sample_ids)
      .def("map_positions", &ViterbiPath::map_positions, py::arg("positions"))
      .def("dump_data_in_range", &ViterbiPath::dump_data_in_range, py::arg("start") = -1,
           py::arg("end") = -1)
      .def("size", &ViterbiPath::size);

  py::class_<MatchGroup>(m, "MatchGroup")
      .def_readonly("match_candidates", &MatchGroup::match_candidates)
      .def_readonly("match_candidates_counts", &MatchGroup::match_candidates_counts)
      .def_readonly("top_four_maps", &MatchGroup::top_four_maps)
      .def_readonly("cm_position", &MatchGroup::cm_position);

  py::class_<Matcher>(m, "Matcher")
      .def(py::init<int, std::vector<double>, double, double, int, int>(), "Initialize",
           py::arg("num_samples"), py::arg("genetic_positions"), py::arg("query_interval_size"),
           py::arg("match_group_interval_size"), py::arg("neighborhood_size") = 4,
           py::arg("min_matches") = 4)
      .def_readonly("query_sites", &Matcher::query_sites)
      .def_readonly("genetic_positions", &Matcher::genetic_positions)
      .def_readonly("query_interval_size", &Matcher::query_interval_size)
      .def_readonly("num_samples", &Matcher::num_samples)
      .def_readonly("num_sites", &Matcher::num_sites)
      .def("process_site", &Matcher::process_site)
      .def("process_all_sites", &Matcher::process_all_sites)
      .def("process_all_sites_numpy", [](Matcher& self, py::array_t<int32_t, py::array::c_style | py::array::forcecast> arr) {
        auto buf = arr.request();
        if (buf.ndim != 2) throw std::runtime_error("Expected 2D array (n_sites × n_haps)");
        int n_sites = static_cast<int>(buf.shape[0]);
        int n_haps = static_cast<int>(buf.shape[1]);
        self.process_all_sites_flat(static_cast<const int32_t*>(buf.ptr), n_sites, n_haps);
      })
      .def("propagate_adjacent_matches", &Matcher::propagate_adjacent_matches)
      .def("get_matches", &Matcher::get_matches)
      .def("serializable_matches", &Matcher::serializable_matches)
      .def("cm_positions", &Matcher::cm_positions)
      .def("get_sorting", &Matcher::get_sorting)
      .def("get_permutation", &Matcher::get_permutation)
      .def("clear", &Matcher::clear);

  py::class_<ImputationMatcher>(m, "ImputationMatcher")
      .def(py::init<int, int, const std::vector<double>&, double, int>(), "Initialize",
           py::arg("num_reference"), py::arg("num_target"), py::arg("genetic_positions"),
           py::arg("query_interval_size"), py::arg("neighborhood_size") = 8)
      .def("process_site", &ImputationMatcher::process_site)
      .def("get_matches", &ImputationMatcher::get_matches);

  py::class_<ThreadsFastLS>(m, "ThreadsFastLS")
      .def(py::init<std::vector<double>, std::vector<double>, double, double, bool, int, bool, int,
                    int>(),
           "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate"), py::arg("Ne") = 2e4, py::arg("sparse_sites") = false,
           py::arg("n_prune") = -1, py::arg("use_hmm") = false, py::arg("burn_in_left") = 0,
           py::arg("burn_in_right") = 0)
      .def(py::init<std::vector<double>, std::vector<double>, double, std::vector<double>,
                    std::vector<double>, bool, int, bool, int, int>(),
           "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate"), py::arg("Ne_sizes"), py::arg("Ne_times"),
           py::arg("sparse_sites") = false, py::arg("n_prune") = -1, py::arg("use_hmm") = false,
           py::arg("burn_in_left") = 0, py::arg("burn_in_right") = 0)
      .def("insert", py::overload_cast<const std::vector<bool>&>(&ThreadsFastLS::insert),
           py::arg("genotypes"))
      .def("impute", &ThreadsFastLS::impute);

  py::class_<ThreadingInstructions>(m, "ThreadingInstructions")
      .def(py::init<const std::vector<ViterbiPath>, const int, const int, const std::vector<int>&>(), "initialize",
           py::arg("paths"), py::arg("start"), py::arg("end"), py::arg("positions"))
      .def(py::init<const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<int>>&,
           const std::vector<std::vector<int>>&, const std::vector<int>&, int, int>(),
           "Initialize", py::arg("starts"), py::arg("tmrcas"), py::arg("targets"), py::arg("mismatches"),
           py::arg("positions"),  py::arg("start"),  py::arg("end"))
      .def_readonly("positions", &ThreadingInstructions::positions)
      .def_readonly("num_sites", &ThreadingInstructions::num_sites)
      .def_readonly("num_samples", &ThreadingInstructions::num_samples)
      .def_readonly("start", &ThreadingInstructions::start)
      .def_readonly("end", &ThreadingInstructions::end)
      .def("all_starts", &ThreadingInstructions::all_starts,
           "Return all segment start positions (bp) as a flat list, concatenated across samples.")
      .def("all_tmrcas", &ThreadingInstructions::all_tmrcas,
           "Return all segment TMRCAs (generations) as a flat list, concatenated across samples.")
      .def("all_targets", &ThreadingInstructions::all_targets,
           "Return all segment target sample indices as a flat list, concatenated across samples.")
      .def("all_mismatches", &ThreadingInstructions::all_mismatches,
           "Return all mismatch site indices as a flat list, concatenated across samples.")
      .def("sub_range", &ThreadingInstructions::sub_range,
           "Return a new ThreadingInstructions restricted to sites in [start_bp, end_bp].")
      .def("add_variants", [](const ThreadingInstructions& self,
              const std::vector<int>& new_positions,
              py::array_t<int, py::array::c_style | py::array::forcecast> geno_arr) {
          auto buf = geno_arr.request();
          if (buf.ndim != 2)
              throw std::runtime_error("new_genotypes must be 2D (n_new_sites × num_samples)");
          int n_new = static_cast<int>(buf.shape[0]);
          int n_samp = static_cast<int>(buf.shape[1]);
          if (n_samp != self.num_samples)
              throw std::runtime_error("new_genotypes columns must equal num_samples");
          if (n_new != static_cast<int>(new_positions.size()))
              throw std::runtime_error("new_positions length must match new_genotypes rows");
          const int* ptr = static_cast<const int*>(buf.ptr);
          std::vector<int> geno_flat(ptr, ptr + static_cast<size_t>(n_new) * n_samp);
          return self.add_variants(new_positions, geno_flat, n_new);
      }, py::arg("new_positions"), py::arg("new_genotypes"))
      .def(py::pickle(
           &threading_instructions_get_state,
           &threading_instructions_set_state))
      .def("left_multiply", &ThreadingInstructions::left_multiply, py::arg("x"), py::arg("diploid") = false, py::arg("normalize") = false,
           py::call_guard<py::gil_scoped_release>(),
           "Compute G^T x (dense). x has length num_samples (or num_samples/2 if diploid). Returns vector of length num_sites.")
      .def("right_multiply", &ThreadingInstructions::right_multiply, py::arg("x"), py::arg("diploid") = false, py::arg("normalize") = false,
           py::call_guard<py::gil_scoped_release>(),
           "Compute G x (dense). x has length num_sites. Returns vector of length num_samples (or num_samples/2 if diploid).")
      .def("materialize_genotypes", &ThreadingInstructions::materialize_genotypes,
           "Cache the dense genotype matrix in memory. Called automatically by multiply methods.")
      .def("genotype_matrix_numpy", [](ThreadingInstructions& self) {
        const auto& gmat = self.get_genotype_matrix();
        int m = self.num_sites;
        int n = self.num_samples;
        py::array_t<int8_t> out({m, n});
        auto* dst = out.mutable_data();
        const int* src = gmat.data();
        size_t total = static_cast<size_t>(m) * n;
        for (size_t i = 0; i < total; i++) dst[i] = static_cast<int8_t>(src[i]);
        return out;
      }, "Return genotype matrix as (n_sites, n_samples) int8 numpy array.")
      .def("materialize_normalized_haploid", &ThreadingInstructions::materialize_normalized_haploid,
           "Cache the mean-centered, variance-normalized haploid genotype matrix.")
      .def("materialize_normalized_diploid", &ThreadingInstructions::materialize_normalized_diploid,
           "Cache the mean-centered, variance-normalized diploid dosage matrix.")
      .def("simplify_for_multiply", &ThreadingInstructions::simplify_for_multiply,
           "Strip TMRCAs and merge consecutive segments with the same target. Reduces multiply cost without changing genotypes.")
      .def("coarsen", &ThreadingInstructions::coarsen, py::arg("min_sites"),
           "Merge segments shorter than min_sites (lossy). Reduces segment count for faster multiply at the cost of genotype accuracy.")
      .def("prepare_rle_multiply", &ThreadingInstructions::prepare_rle_multiply,
           "Precompute run-length encoding of 1-runs per sample. Sparse baseline multiply. "
           "O(m + total_1_runs) per multiply; prefer DAG multiply for large n.")
      .def("right_multiply_rle", [](ThreadingInstructions& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> arr,
              bool diploid, bool normalize) {
        auto buf = arr.request();
        const int m = self.num_sites;
        const int n = self.num_samples;
        NormStats ns;
        if (normalize) ns = get_norm_stats(self, diploid);
        int k;
        std::vector<double> result;
        if (buf.ndim == 1) {
            k = 1;
            if (static_cast<int>(buf.shape[0]) != m)
                throw std::runtime_error("1D input must have length num_sites");
            std::vector<double> x(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + m);
            std::vector<double> offsets;
            if (normalize) norm_prescale_right(x, ns, offsets, m, 1);
            result = self.right_multiply_rle(x);
            if (normalize) norm_postprocess_right(result, offsets, diploid, n, 1);
            else if (diploid) diploid_reduce_right(result, n, 1);
        } else if (buf.ndim == 2) {
            if (static_cast<int>(buf.shape[0]) != m)
                throw std::runtime_error("First dimension must be num_sites");
            k = static_cast<int>(buf.shape[1]);
            std::vector<double> x(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + buf.size);
            std::vector<double> offsets;
            if (normalize) norm_prescale_right(x, ns, offsets, m, k);
            result = self.right_multiply_rle_batch(x, k);
            if (normalize) norm_postprocess_right(result, offsets, diploid, n, k);
            else if (diploid) diploid_reduce_right(result, n, k);
        } else {
            throw std::runtime_error("Input must be 1D or 2D");
        }
        int out_rows = (diploid ? n / 2 : n);
        if (k == 1) {
            py::array_t<double> out(static_cast<size_t>(out_rows));
            std::memcpy(out.mutable_data(), result.data(), out_rows * sizeof(double));
            return out;
        }
        py::array_t<double> out({out_rows, k});
        std::memcpy(out.mutable_data(), result.data(), result.size() * sizeof(double));
        return out;
      }, py::arg("X"), py::arg("diploid")=false, py::arg("normalize")=false,
           "Compute G @ X via RLE. Supports diploid and normalize flags.")
      .def("left_multiply_rle", [](ThreadingInstructions& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> arr,
              bool diploid, bool normalize) {
        auto buf = arr.request();
        const int m = self.num_sites;
        const int n = self.num_samples;
        const int n_dip = n / 2;
        const int expected_rows = diploid ? n_dip : n;
        NormStats ns;
        if (normalize) ns = get_norm_stats(self, diploid);
        int k;
        std::vector<double> result;
        if (buf.ndim == 1) {
            k = 1;
            if (static_cast<int>(buf.shape[0]) != expected_rows)
                throw std::runtime_error("1D input length mismatch");
            std::vector<double> y(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + expected_rows);
            auto sums = col_sums_of(y.data(), expected_rows, 1);
            if (diploid) y = diploid_expand_left(y, n, 1);
            result = self.left_multiply_rle(y);
            if (normalize) norm_postprocess_left(result, ns, sums, m, 1);
        } else if (buf.ndim == 2) {
            if (static_cast<int>(buf.shape[0]) != expected_rows)
                throw std::runtime_error("First dimension mismatch");
            k = static_cast<int>(buf.shape[1]);
            std::vector<double> y(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + buf.size);
            auto sums = col_sums_of(y.data(), expected_rows, k);
            if (diploid) y = diploid_expand_left(y, n, k);
            result = self.left_multiply_rle_batch(y, k);
            if (normalize) norm_postprocess_left(result, ns, sums, m, k);
        } else {
            throw std::runtime_error("Input must be 1D or 2D");
        }
        if (k == 1) {
            py::array_t<double> out(static_cast<size_t>(m));
            std::memcpy(out.mutable_data(), result.data(), m * sizeof(double));
            return out;
        }
        py::array_t<double> out({m, k});
        std::memcpy(out.mutable_data(), result.data(), result.size() * sizeof(double));
        return out;
      }, py::arg("X"), py::arg("diploid")=false, py::arg("normalize")=false,
           "Compute G.T @ X via RLE. Supports diploid and normalize flags.")
      .def("het_per_site", &ThreadingInstructions::het_per_site,
           "Count heterozygous diploid individuals per site. Returns int32 array "
           "of length num_sites. O(n*m) time, O(n) memory (no dense matrix allocated).")
      .def("het_per_individual", &ThreadingInstructions::het_per_individual,
           "Count heterozygous sites per diploid individual. Returns int32 array "
           "of length num_samples/2. O(n*m) time, O(n) memory.")
      .def_property_readonly("dag_num_sets", &ThreadingInstructions::get_dag_num_sets,
           "Total mutation set nodes across all DAG chunks (0 until prepare_dag_multiply is called).")
      .def_property_readonly("dag_total_edges", &ThreadingInstructions::get_dag_total_edges)
      .def_property_readonly("dag_total_connections", &ThreadingInstructions::get_dag_total_connections)
      .def_property_readonly("dag_total_mutations", &ThreadingInstructions::get_dag_total_mutations)
      .def_property_readonly("dag_num_chunks", &ThreadingInstructions::get_dag_num_chunks)
      .def("prepare_dag_multiply", &ThreadingInstructions::prepare_dag_multiply,
           py::call_guard<py::gil_scoped_release>(),
           py::arg("num_chunks") = 0,
           "Precompute region-chunked mutation-set DAGs for sublinear multiply, "
           "following the ARG multiplication approach used in arg-needle-lib. "
           "num_chunks=0 (default) uses OMP_NUM_THREADS. Each chunk covers a disjoint "
           "site range and is processed in parallel during multiply.")
      .def("right_multiply_dag", [](ThreadingInstructions& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> arr,
              bool diploid, bool normalize) {
        auto buf = arr.request();
        const int m = self.num_sites;
        const int n = self.num_samples;
        NormStats ns;
        if (normalize) ns = get_norm_stats(self, diploid);
        int k;
        std::vector<double> result;
        if (buf.ndim == 1) {
            k = 1;
            if (static_cast<int>(buf.shape[0]) != m)
                throw std::runtime_error("1D input must have length num_sites");
            std::vector<double> x(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + m);
            std::vector<double> offsets;
            if (normalize) norm_prescale_right(x, ns, offsets, m, 1);
            { py::gil_scoped_release release; result = self.right_multiply_dag(x); }
            if (normalize) norm_postprocess_right(result, offsets, diploid, n, 1);
            else if (diploid) diploid_reduce_right(result, n, 1);
        } else if (buf.ndim == 2) {
            if (static_cast<int>(buf.shape[0]) != m)
                throw std::runtime_error("First dimension must be num_sites");
            k = static_cast<int>(buf.shape[1]);
            std::vector<double> x(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + buf.size);
            std::vector<double> offsets;
            if (normalize) norm_prescale_right(x, ns, offsets, m, k);
            { py::gil_scoped_release release; result = self.right_multiply_dag_batch(x, k); }
            if (normalize) norm_postprocess_right(result, offsets, diploid, n, k);
            else if (diploid) diploid_reduce_right(result, n, k);
        } else {
            throw std::runtime_error("Input must be 1D or 2D");
        }
        int out_rows = (diploid ? n / 2 : n);
        if (k == 1) {
            py::array_t<double> out(static_cast<size_t>(out_rows));
            std::memcpy(out.mutable_data(), result.data(), out_rows * sizeof(double));
            return out;
        }
        py::array_t<double> out({out_rows, k});
        std::memcpy(out.mutable_data(), result.data(), result.size() * sizeof(double));
        return out;
      }, py::arg("X"), py::arg("diploid")=false, py::arg("normalize")=false,
           "Compute G @ X via mutation-set DAG (arg-needle-lib approach). Supports diploid and normalize flags. Sublinear in n.")
      .def("left_multiply_dag", [](ThreadingInstructions& self,
              py::array_t<double, py::array::c_style | py::array::forcecast> arr,
              bool diploid, bool normalize) {
        auto buf = arr.request();
        const int m = self.num_sites;
        const int n = self.num_samples;
        const int n_dip = n / 2;
        const int expected_rows = diploid ? n_dip : n;
        NormStats ns;
        if (normalize) ns = get_norm_stats(self, diploid);
        int k;
        std::vector<double> result;
        if (buf.ndim == 1) {
            k = 1;
            if (static_cast<int>(buf.shape[0]) != expected_rows)
                throw std::runtime_error("1D input length mismatch");
            std::vector<double> y(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + expected_rows);
            auto sums = col_sums_of(y.data(), expected_rows, 1);
            if (diploid) y = diploid_expand_left(y, n, 1);
            { py::gil_scoped_release release; result = self.left_multiply_dag(y); }
            if (normalize) norm_postprocess_left(result, ns, sums, m, 1);
        } else if (buf.ndim == 2) {
            if (static_cast<int>(buf.shape[0]) != expected_rows)
                throw std::runtime_error("First dimension mismatch");
            k = static_cast<int>(buf.shape[1]);
            std::vector<double> y(static_cast<double*>(buf.ptr), static_cast<double*>(buf.ptr) + buf.size);
            auto sums = col_sums_of(y.data(), expected_rows, k);
            if (diploid) y = diploid_expand_left(y, n, k);
            { py::gil_scoped_release release; result = self.left_multiply_dag_batch(y, k); }
            if (normalize) norm_postprocess_left(result, ns, sums, m, k);
        } else {
            throw std::runtime_error("Input must be 1D or 2D");
        }
        if (k == 1) {
            py::array_t<double> out(static_cast<size_t>(m));
            std::memcpy(out.mutable_data(), result.data(), m * sizeof(double));
            return out;
        }
        py::array_t<double> out({m, k});
        std::memcpy(out.mutable_data(), result.data(), result.size() * sizeof(double));
        return out;
      }, py::arg("X"), py::arg("diploid")=false, py::arg("normalize")=false,
           "Compute G.T @ X via mutation-set DAG (arg-needle-lib approach). Supports diploid and normalize flags. Sublinear in n.")
      .def("visit_clades", [](const ThreadingInstructions& self, py::function callback) {
        self.visit_clades([&](const std::vector<int>& descendants, double height,
                              int start_bp, int end_bp) {
          py::gil_scoped_acquire acquire;
          callback(descendants, height, start_bp, end_bp);
        });
      }, py::arg("callback"),
           "Traverse all clades in the implicit ARG. Calls callback(descendants, height, start_bp, end_bp) for each clade.")
      .def("visit_branches", [](const ThreadingInstructions& self, py::function callback) {
        self.visit_branches([&](const std::vector<int>& descendants, double child_height,
                                double parent_height, int start_bp, int end_bp) {
          py::gil_scoped_acquire acquire;
          callback(descendants, child_height, parent_height, start_bp, end_bp);
        });
      }, py::arg("callback"),
           "Traverse all branches in the implicit ARG. Calls callback(descendants, child_height, parent_height, start_bp, end_bp).")
      .def("total_volume", &ThreadingInstructions::total_volume,
           "Total branch volume (sum of branch_length * span over all branches). Uses lightweight union-find (no member lists).")
      .def("allele_frequency_spectrum", &ThreadingInstructions::allele_frequency_spectrum,
           "Branch-length AFS: afs[k] = total volume of branches with k descendants. Returns list of length num_samples+1. "
           "Uses count-only union-find; ~50-70x faster than arg-needle-lib's bitset_volume_map.")
      .def("association_diploid", [](const ThreadingInstructions& self,
              const std::vector<double>& phenotypes, double p_threshold) {
        auto r = self.association_diploid(phenotypes, p_threshold);
        return py::make_tuple(r.position, r.chi2, r.pvalue, r.n_descendants);
      }, py::arg("phenotypes"), py::arg("p_threshold") = 5e-8,
           "Clade-based association testing. Tests each clade's diploid dosage against phenotypes via chi-square. "
           "Returns (positions, chi2_stats, p_values, n_descendants) for clades passing p_threshold.")
      .def("generate_mutations", [](const ThreadingInstructions& self,
              double mutation_rate, int seed) {
        auto r = self.generate_mutations(mutation_rate, seed);
        int n = self.num_samples;
        int m = r.n_mutations;
        py::array_t<int> positions(m);
        std::memcpy(positions.mutable_data(), r.positions.data(), m * sizeof(int));
        py::array_t<int> genotypes({m, n});
        std::memcpy(genotypes.mutable_data(), r.genotypes.data(), m * n * sizeof(int));
        return py::make_tuple(positions, genotypes);
      }, py::arg("mutation_rate"), py::arg("seed") = 42,
           "Poisson-sample mutations on branches proportional to volume. "
           "Returns (positions_array, genotypes_matrix) where genotypes is (n_mutations, num_samples) int32.");

  py::class_<ConsistencyWrapper>(m, "ConsistencyWrapper")
      .def(py::init<const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<int>>&,
           const std::vector<std::vector<int>>&, const std::vector<int>&, const std::vector<double>&>(),
           "Initialize", py::arg("starts"), py::arg("tmrcas"), py::arg("targets"), py::arg("mismatches"),
           py::arg("physical_positions"), py::arg("allele_ages"))
      .def(py::init<ThreadingInstructions&, const std::vector<double>&>(),
           "Initialize", py::arg("instructions"), py::arg("allele_ages"))
      .def("process_site", &ConsistencyWrapper::process_site)
      .def("get_consistent_instructions", &ConsistencyWrapper::get_consistent_instructions);

  m.def("run_consistency", &run_consistency, py::arg("instructions"), py::arg("allele_ages"));

  py::class_<AgeEstimator>(m, "AgeEstimator")
      .def(py::init<const ThreadingInstructions&>(), "initialize", py::arg("instructions"))
      .def("process_site", &AgeEstimator::process_site)
      .def("get_inferred_ages", &AgeEstimator::get_inferred_ages);

  m.def("estimate_ages", &estimate_ages, py::arg("instructions"));

  py::class_<GenotypeIterator>(m, "GenotypeIterator")
      .def(py::init<const ThreadingInstructions&>(), "initialize", py::arg("instructions"))
      .def("next_genotype", &GenotypeIterator::next_genotype)
      .def("has_next_genotype", &GenotypeIterator::has_next_genotype);

  py::class_<VCFWriter>(m, "VCFWriter")
      .def(py::init<ThreadingInstructions&>(), "initialize", py::arg("instructions"))
      .def("set_chrom", &VCFWriter::set_chrom)
      .def("set_pos", &VCFWriter::set_pos)
      .def("set_id", &VCFWriter::set_id)
      .def("set_ref", &VCFWriter::set_ref)
      .def("set_alt", &VCFWriter::set_alt)
      .def("set_qual", &VCFWriter::set_qual)
      .def("set_filter", &VCFWriter::set_filter)
      .def("set_sample_names", &VCFWriter::set_sample_names)
      .def("write_vcf", &VCFWriter::write_vcf);

  // Forward-backward Li-Stephens algorithm (replaces numba JIT'd fwbw)
  m.def("forwards_ls_hap", [](
      py::array_t<double, py::array::c_style | py::array::forcecast> H,
      py::array_t<double, py::array::c_style | py::array::forcecast> s,
      py::array_t<double, py::array::c_style | py::array::forcecast> e,
      py::array_t<double, py::array::c_style | py::array::forcecast> r) {
    auto H_buf = H.request();
    auto s_buf = s.request();
    int m_sites = H_buf.shape[0];
    int n_refs = H_buf.shape[1];

    auto [F_vec, c_vec] = forwards_ls_hap(
        n_refs, m_sites,
        static_cast<const double*>(H_buf.ptr),
        static_cast<const double*>(s_buf.ptr),
        static_cast<const double*>(e.request().ptr),
        static_cast<const double*>(r.request().ptr));

    py::array_t<double> F({m_sites, n_refs});
    std::memcpy(F.mutable_data(), F_vec.data(), F_vec.size() * sizeof(double));
    py::array_t<double> c(m_sites);
    std::memcpy(c.mutable_data(), c_vec.data(), c_vec.size() * sizeof(double));
    return py::make_tuple(F, c);
  });

  m.def("backwards_ls_hap", [](
      py::array_t<double, py::array::c_style | py::array::forcecast> H,
      py::array_t<double, py::array::c_style | py::array::forcecast> s,
      py::array_t<double, py::array::c_style | py::array::forcecast> e,
      py::array_t<double, py::array::c_style | py::array::forcecast> c,
      py::array_t<double, py::array::c_style | py::array::forcecast> r) {
    auto H_buf = H.request();
    int m_sites = H_buf.shape[0];
    int n_refs = H_buf.shape[1];

    auto B_vec = backwards_ls_hap(
        n_refs, m_sites,
        static_cast<const double*>(H_buf.ptr),
        static_cast<const double*>(s.request().ptr),
        static_cast<const double*>(e.request().ptr),
        static_cast<const double*>(c.request().ptr),
        static_cast<const double*>(r.request().ptr));

    py::array_t<double> B({m_sites, n_refs});
    std::memcpy(B.mutable_data(), B_vec.data(), B_vec.size() * sizeof(double));
    return B;
  });

  // .threads file I/O (replaces h5py)
  m.def("serialize_threads", &serialize_threads,
        py::arg("filename"), py::arg("instructions"),
        py::arg("metadata_cols") = std::vector<std::vector<std::string>>(),
        py::arg("allele_ages") = std::vector<double>(),
        py::arg("sample_names") = std::vector<std::string>());
  m.def("deserialize_threads", &deserialize_threads, py::arg("filename"));
  m.def("read_threads_metadata", &read_threads_metadata, py::arg("filename"));
  m.def("read_threads_sample_names", &read_threads_sample_names, py::arg("filename"));
}
