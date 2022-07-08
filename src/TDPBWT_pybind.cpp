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
        .def_readonly("site", &Node::site)
        .def_readonly("genotype", &Node::genotype);

    // py::class_<TDPBWT>(m, "TDPBWT")
    //     .def(py::init<std::vector<int>, std::vector<float>, float>(), "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
    //        py::arg("mutation_rate") = 0)
    //     .def_readonly("num_samples", &TDPBWT::num_samples)
    //     .def_readonly("num_sites", &TDPBWT::num_sites)
    //     .def_readonly("mutation_rate", &TDPBWT::mutation_rate)
    //     .def("insert", &TDPBWT::insert)
    //     .def("longest_suffix", &TDPBWT::longest_suffix)

    py::class_<DPBWT>(m, "DPBWT")
        .def(py::init<std::vector<double>, std::vector<double>, double>(), "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate") = 0)
        .def_readonly("num_samples", &DPBWT::num_samples)
        .def_readonly("num_sites", &DPBWT::num_sites)
        .def_readonly("mutation_rate", &DPBWT::mutation_rate)
        .def("insert", &DPBWT::insert)
        .def("print_sorting", &DPBWT::print_sorting)
        .def("print_divergence", &DPBWT::print_divergence)
        .def("longest_prefix", &DPBWT::longest_prefix)
        .def("fastLS", &DPBWT::fastLS);
}
