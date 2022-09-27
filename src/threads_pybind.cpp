#include "DPBWT.hpp"

// #include <deque>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <sstream>
// #include <stdexcept>
#include <vector>

namespace py = pybind11;
// using std::deque;
// using std::vector;

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

    py::class_<DPBWT>(m, "DPBWT")
        .def(py::init<std::vector<double>, std::vector<double>, double, double, bool, int, bool>(), "Initialize",
             py::arg("physical_positions"), py::arg("genetic_positions"),
             py::arg("mutation_rate"), py::arg("Ne") = 2e4, py::arg("sparse_sites")=false, py::arg("n_prune")=-1, py::arg("use_hmm")=false)
        .def(py::init<std::vector<double>, std::vector<double>, double, std::vector<double>, std::vector<double>, bool, int, bool>(), "Initialize",
             py::arg("physical_positions"), py::arg("genetic_positions"),
             py::arg("mutation_rate"), py::arg("Ne_sizes"), py::arg("Ne_times"),
             py::arg("sparse_sites")=false, py::arg("n_prune")=-1, py::arg("use_hmm")=false)
        .def_readonly("num_samples", &DPBWT::num_samples)
        .def_readonly("num_sites", &DPBWT::num_sites)
        .def_readonly("mutation_rate", &DPBWT::mutation_rate)
        .def_readonly("bp_sizes", &DPBWT::bp_sizes)
        .def_readonly("cm_sizes", &DPBWT::cm_sizes)
        .def_readonly("mutation_rate", &DPBWT::mutation_rate)
        .def_readonly("n_prune", &DPBWT::mutation_rate)
        .def_readonly("ID_map", &DPBWT::ID_map)
        .def_readonly("demography", &DPBWT::demography)
        .def_readonly("hmm", &DPBWT::hmm)
        .def("delete_hmm", &DPBWT::delete_hmm)
        .def("mutation_penalties", &DPBWT::mutation_penalties)
        .def("recombination_penalties", &DPBWT::recombination_penalties)
        .def("date_segment", &DPBWT::date_segment, py::arg("id1"), py::arg("id2"), py::arg("start"), py::arg("end"))
        .def("insert", py::overload_cast<const std::vector<bool>&>(&DPBWT::insert), py::arg("genotypes"))
        .def("insert", py::overload_cast<const int, const std::vector<bool>&>(&DPBWT::insert), py::arg("ID"), py::arg("genotypes"))
        .def("delete_ID", &DPBWT::delete_ID, py::arg("ID"))
        .def("print_sorting", &DPBWT::print_sorting)
        .def("thread", py::overload_cast<const std::vector<bool>&>(&DPBWT::thread), py::arg("genotypes"))
        .def("thread", py::overload_cast<const int, const std::vector<bool>&>(&DPBWT::thread), py::arg("new_sample_ID"), py::arg("genotypes"))
        .def("thread_with_mutations", py::overload_cast<const std::vector<bool>&>(&DPBWT::thread_with_mutations), py::arg("genotypes"))
        .def("thread_with_mutations", py::overload_cast<const int, const std::vector<bool>&>(&DPBWT::thread_with_mutations), py::arg("new_sample_ID"), py::arg("genotypes"))
        // .def("thread_from_file", &DPBWT::thread_from_file, py::arg("hack_gz"), py::arg("n_cycle")=0)
        .def("fastLS", &DPBWT::fastLS);
}
