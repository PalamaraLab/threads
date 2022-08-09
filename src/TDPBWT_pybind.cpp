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

PYBIND11_MODULE(TDPBWT_python_bindings, m) {
    py::class_<Node>(m, "Node")
        .def_readonly("sample_ID", &Node::sample_ID)
        .def_readonly("genotype", &Node::genotype);
    py::class_<Demography>(m, "Demography")
        .def_readonly("times", &Demography::times)
        .def_readonly("sizes", &Demography::sizes)
        .def("std_to_gen", &Demography::std_to_gen);

    // py::class_<TDPBWT>(m, "TDPBWT")
    //     .def(py::init<std::vector<int>, std::vector<float>, float>(), "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
    //        py::arg("mutation_rate") = 0)
    //     .def_readonly("num_samples", &TDPBWT::num_samples)
    //     .def_readonly("num_sites", &TDPBWT::num_sites)
    //     .def_readonly("mutation_rate", &TDPBWT::mutation_rate)
    //     .def("insert", &TDPBWT::insert)
    //     .def("longest_suffix", &TDPBWT::longest_suffix)

    py::class_<DPBWT>(m, "DPBWT")
        .def(py::init<std::vector<double>, std::vector<double>, double, double, std::string>(), "Initialize",
        py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate"), py::arg("Ne") = 2e4, py::arg("mode")="ML")
        .def(py::init<std::vector<double>, std::vector<double>, double, std::vector<double>, std::vector<double>, std::string>(), "Initialize",
        py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate"), py::arg("Ne_sizes"), py::arg("Ne_times"), py::arg("mode")="ML")
        .def_readonly("num_samples", &DPBWT::num_samples)
        .def_readonly("num_sites", &DPBWT::num_sites)
        .def_readonly("mutation_rate", &DPBWT::mutation_rate)
        .def_readonly("bp_sizes", &DPBWT::bp_sizes)
        .def_readonly("cm_sizes", &DPBWT::cm_sizes)
        .def_readonly("mutation_rate", &DPBWT::mutation_rate)
        .def_readonly("ID_map", &DPBWT::ID_map)
        .def_readonly("demography", &DPBWT::demography)
        .def("mutation_penalties", &DPBWT::mutation_penalties)
        .def("recombination_penalties", &DPBWT::recombination_penalties)
        .def("date_segment", &DPBWT::date_segment, py::arg("id1"), py::arg("id2"), py::arg("start"), py::arg("end"))
        .def("insert", py::overload_cast<std::vector<bool>>(&DPBWT::insert), py::arg("genotypes"))
        .def("insert", py::overload_cast<int, std::vector<bool>>(&DPBWT::insert), py::arg("ID"), py::arg("genotypes"))
        .def("delete_ID", &DPBWT::delete_ID, py::arg("ID"))
        .def("print_sorting", &DPBWT::print_sorting)
        .def("thread", py::overload_cast<std::vector<bool>>(&DPBWT::thread), py::arg("genotypes"))
        .def("thread", py::overload_cast<int, std::vector<bool>>(&DPBWT::thread), py::arg("new_sample_ID"), py::arg("genotypes"))
        .def("fastLS", &DPBWT::fastLS);
}
