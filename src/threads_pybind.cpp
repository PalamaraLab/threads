#include "TGEN.hpp"
#include "Threads.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;
PYBIND11_MODULE(threads_python_bindings, m) {
  py::class_<Node>(m, "Node")
      .def_readonly("sample_ID", &Node::sample_ID)
      .def_readonly("genotype", &Node::genotype);
  py::class_<Demography>(m, "Demography")
      .def_readonly("times", &Demography::times)
      .def_readonly("sizes", &Demography::sizes)
      .def("std_to_gen", &Demography::std_to_gen)
      .def("expected_branch_length", &Demography::expected_branch_length);

  py::class_<HMM>(m, "HMM")
      .def_readonly("num_states", &HMM::num_states)
      .def_readonly("expected_times", &HMM::expected_times)
      .def_readonly("trellis", &HMM::trellis)
      .def_readonly("pointers", &HMM::pointers)
      .def_readonly("hom_score", &HMM::hom_score)
      .def_readonly("het_score", &HMM::het_score)
      .def_readonly("transition_score", &HMM::transition_score)
      .def_readonly("non_transition_score", &HMM::non_transition_score)
      .def("breakpoints", &HMM::breakpoints);

  py::class_<Threads>(m, "Threads")
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
      .def_readonly("num_samples", &Threads::num_samples)
      .def_readonly("num_sites", &Threads::num_sites)
      .def_readonly("mutation_rate", &Threads::mutation_rate)
      .def_readonly("bp_sizes", &Threads::bp_sizes)
      .def_readonly("cm_sizes", &Threads::cm_sizes)
      .def_readonly("mutation_rate", &Threads::mutation_rate)
      .def_readonly("n_prune", &Threads::mutation_rate)
      .def_readonly("ID_map", &Threads::ID_map)
      .def_readonly("demography", &Threads::demography)
      .def_readonly("hmm", &Threads::hmm)
      .def_readonly("threading_start", &Threads::threading_start)
      .def_readonly("threading_end", &Threads::threading_end)
      .def("trimmed_positions", &Threads::trimmed_positions)
      .def("set_impute", &Threads::set_impute, py::arg("impute_state"))
      .def("delete_hmm", &Threads::delete_hmm)
      .def("mutation_penalties", &Threads::mutation_penalties)
      .def("recombination_penalties", &Threads::recombination_penalties)
      .def("date_segment", &Threads::date_segment, py::arg("id1"), py::arg("id2"), py::arg("start"),
           py::arg("end"))
      .def("insert", py::overload_cast<const std::vector<bool>&>(&Threads::insert),
           py::arg("genotypes"))
      .def("insert", py::overload_cast<const int, const std::vector<bool>&>(&Threads::insert),
           py::arg("ID"), py::arg("genotypes"))
      .def("remove", &Threads::remove, py::arg("ID"))
      .def("print_sorting", &Threads::print_sorting)
      .def("thread", py::overload_cast<const std::vector<bool>&, const int>(&Threads::thread),
           py::arg("genotypes"), py::arg("L") = 1)
      .def("thread",
           py::overload_cast<const int, const std::vector<bool>&, const int>(&Threads::thread),
           py::arg("new_sample_ID"), py::arg("genotypes"), py::arg("L") = 1)
      // .def("thread_with_mutations", py::overload_cast<const
      // std::vector<bool>&>(&Threads::thread_with_mutations), py::arg("genotypes"))
      // .def("thread_with_mutations", py::overload_cast<const int, const
      // std::vector<bool>&>(&Threads::thread_with_mutations), py::arg("new_sample_ID"),
      // py::arg("genotypes")) .def("thread_from_file", &Threads::thread_from_file,
      // py::arg("hack_gz"), py::arg("n_cycle")=0)
      .def("fastLS", &Threads::fastLS);

  py::class_<TGEN>(m, "TGEN")
      .def(py::init<std::vector<int>, std::vector<std::vector<int>>, std::vector<std::vector<int>>,
                    std::vector<std::vector<int>>>(),
           "Initialize", py::arg("positions"), py::arg("bp_starts"), py::arg("target_ids"),
           py::arg("het_sites"))
      .def("query", &TGEN::query, py::return_value_policy::reference_internal);
}
