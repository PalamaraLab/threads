// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;
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
      .def("process_site_hets", &ThreadsLowMem::process_site_hets)
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
      .def("all_starts", &ThreadingInstructions::all_starts)
      .def("all_tmrcas", &ThreadingInstructions::all_tmrcas)
      .def("all_targets", &ThreadingInstructions::all_targets)
      .def("all_mismatches", &ThreadingInstructions::all_mismatches);

  py::class_<ConsistencyWrapper>(m, "ConsistencyWrapper")
      .def(py::init<const std::vector<std::vector<int>>&, const std::vector<std::vector<double>>&, const std::vector<std::vector<int>>&,
           const std::vector<std::vector<int>>&, const std::vector<int>&, const std::vector<double>&>(),
           "Initialize", py::arg("starts"), py::arg("tmrcas"), py::arg("targets"), py::arg("mismatches"),
           py::arg("physical_positions"), py::arg("allele_ages"))
      .def(py::init<ThreadingInstructions&, const std::vector<double>&>(),
           "Initialize", py::arg("instructions"), py::arg("allele_ages"))
      .def("process_site", &ConsistencyWrapper::process_site)
      .def("get_consistent_instructions", &ConsistencyWrapper::get_consistent_instructions);

  py::class_<AgeEstimator>(m, "AgeEstimator")
      .def(py::init<const ThreadingInstructions&>(), "initialize", py::arg("instructions"))
      .def("process_site", &AgeEstimator::process_site)
      .def("get_inferred_ages", &AgeEstimator::get_inferred_ages);

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
}
