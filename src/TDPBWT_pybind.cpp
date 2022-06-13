#include "TDPBWT.hpp"

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

    py::class_<TDPBWT>(m, "TDPBWT")
        .def(py::init<std::vector<int>, std::vector<float>, float>(), "Initialize", py::arg("physical_positions"), py::arg("genetic_positions"),
           py::arg("mutation_rate") = 0)
        .def_readonly("num_samples", &TDPBWT::num_samples)
        .def_readonly("num_sites", &TDPBWT::num_sites)
        .def_readonly("mutation_rate", &TDPBWT::mutation_rate)
        .def("insert", &TDPBWT::insert)
        .def("longest_suffix", &TDPBWT::longest_suffix)
        .def("print_sorting", &TDPBWT::print_sorting);
}

//   int sample_ID;
//   int site;
//   bool genotype;

// PYBIND11_MODULE(arg_needle_python_bindings, m) {
//   py::class_<ARGNode>(m, "ARGNode")
//       .def(py::init<int, real, real, real>(), py::arg("ID"), py::arg("height"), py::arg("start"),
//            py::arg("end"))
//       .def_readonly("ID", &ARGNode::ID)
//       .def_readonly("height", &ARGNode::height)
//       .def_readonly("start", &ARGNode::start)
//       .def_readonly("end", &ARGNode::end)
//       .def("add_parent", &ARGNode::add_parent, py::arg("start"), py::arg("end"), py::arg("parent"))
//       .def("remove_parent", &ARGNode::remove_parent, py::arg("start"))
//       .def("update_parent_start", &ARGNode::update_parent_start, py::arg("start_old"),
//            py::arg("start_new"))
//       .def("update_parent_end", &ARGNode::update_parent_end, py::arg("start"), py::arg("end_new"))
//       // .def_readwrite("parents", &ARGNode::parents)
//       // .def_readonly("parents", &ARGNode::parents, py::return_value_policy::copy)
//       .def(
//           "parent_edges",
//           [](const ARGNode& node) {
//             vector<ARGEdge> edges;
//             for (auto const& map_entry : node.parents) {
//               edges.emplace_back(*(map_entry.second));
//             }
//             return edges;
//           },
//           "Returns a sorted list of parent edges", py::return_value_policy::reference)
//       // .def("parent_edge_at_start",
//       //   [](const ARGNode &node, real start) {
//       //     // TODO: custom exception if not found
//       //     return *node.parents.at(start);
//       //   }, py::return_value_policy::reference)
//       .def("parent_edge_at", &ARGNode::parent_edge_at, py::arg("position"),
//            py::return_value_policy::reference)
//       .def(
//           "parent_starts",
//           [](const ARGNode& node) {
//             vector<real> starts;
//             for (auto const& map_entry : node.parents) {
//               starts.push_back(map_entry.first);
//             }
//             return starts;
//           },
//           "Returns a sorted list of parent starts")
//       // This function commented out by Arni 22/04/22: copy construction of
//       // edges is no longer possible with unique_pointer mutation vectors
//       // .def(
//       //     "child_edges_at", // creates a copy
//       //     [](ARGNode& node, real position) {
//       //       vector<ARGEdge> edges;
//       //       for (auto edge : node.children_at(position)) {
//       //         edges.push_back(*edge);
//       //       }
//       //       return edges;
//       //     },
//       //     py::arg("position"))
//       // This function commented out by Arni 22/04/22: copy construction of
//       // edges is no longer possible with unique_pointer mutation vectors
//       // .def(
//       //     "child_edges_overlap", // creates a copy
//       //     [](ARGNode& node, real pos_start, real pos_end) {
//       //       vector<ARGEdge> edges;
//       //       for (auto edge : node.children_overlap(pos_start, pos_end)) {
//       //         edges.push_back(*edge);
//       //       }
//       //       return edges;
//       //     },
//       //     py::arg("pos_start"), py::arg("pos_end"))
//       .def("__repr__", [](const ARGNode& node) {
//         std::ostringstream oss;
//         oss << node;
//         return oss.str();
//       });

//   py::class_<ARGEdge>(m, "ARGEdge")
//       .def_readonly("start", &ARGEdge::start)
//       .def_readonly("end", &ARGEdge::end)
//       .def(
//           "mutations",
//           [](const ARGEdge& edge) {
//             vector<Mutation> mutations;
//             if (edge.mutations != nullptr) {
//               for (Mutation* mutation : *edge.mutations) {
//                 mutations.push_back(*mutation);
//               }
//             }
//             return mutations;
//           },
//           "Return list of all mutation objects")
//       .def_readwrite("child", &ARGEdge::child, py::return_value_policy::reference)
//       .def_readwrite("parent", &ARGEdge::parent, py::return_value_policy::reference)
//       .def("__repr__", [](const ARGEdge& edge) {
//         std::ostringstream oss;
//         oss << edge;
//         return oss.str();
//       });

//   py::class_<Site>(m, "Site")
//       .def_readonly("ID", &Site::ID)
//       .def_readonly("position", &Site::position);

//   py::class_<Mutation>(m, "Mutation")
//       .def_readonly("edge", &Mutation::edge)
//       .def_readonly("position", &Mutation::position)
//       .def_readonly("height", &Mutation::height)
//       .def_readonly("site", &Mutation::site);
//   // .def("__repr__", [](const Mutation& mut) {
//   //   std::ostringstream oss;
//   //   oss << mut;
//   //   return oss.str();
//   // });

//   py::class_<Root>(m, "Root")
//       .def_readonly("start", &Root::start)
//       .def_readonly("end", &Root::end)
//       .def_readwrite("node", &Root::node, py::return_value_policy::reference)
//       .def("__repr__", [](const Root& root) {
//         std::ostringstream oss;
//         oss << "[" << root.start << ", " << root.end << ")" << endl;
//         oss << *root.node << endl;
//         return oss.str();
//       });

//   py::class_<ARG>(m, "ARG")
//       .def(py::init<real, real, int>(), "Construct an empty ARG", py::arg("start"), py::arg("end"),
//            py::arg("reserved_samples") = 0)
//       .def(py::init<real, real, vector<real>, deque<bool>, vector<std::pair<int, int>>,
//                     vector<std::pair<real, real>>, int>(),
//            "Construct an ARG with node / edge information", py::arg("start"), py::arg("end"),
//            py::arg("node_heights"), py::arg("is_sample"), py::arg("edge_ids"),
//            py::arg("edge_ranges"), py::arg("reserved_samples") = -1)
//       .def(py::init<real, real, int, int, int>(), "Construct a test ARG for memory footprint",
//            py::arg("start"), py::arg("end"), py::arg("num_samples"), py::arg("num_trees"),
//            py::arg("reserved_samples") = -1)
//       .def_readonly("start", &ARG::start)
//       .def_readonly("end", &ARG::end)
//       .def_readonly("threaded_samples", &ARG::threaded_samples)
//       .def_readonly("offset", &ARG::offset)                     // set using set_offset
//       .def_readonly("chromosome", &ARG::chromosome)             // set using set_chromosome
//       .def_readonly("reserved_samples", &ARG::reserved_samples) // set using constructors
//       .def(
//           "node",
//           [](const ARG& arg, int ID) {
//             // TODO: custom exception if not found
//             return arg.arg_nodes.at(ID).get();
//           },
//           py::return_value_policy::reference, "Get node at ID", py::arg("ID"))
//       .def(
//           "node_ids", // unsorted
//           [](const ARG& arg) {
//             vector<int> ids;
//             for (auto const& map_entry : arg.arg_nodes) {
//               ids.push_back(map_entry.first);
//             }
//             return ids;
//           },
//           "Return list of all node IDs")
//       .def(
//           "mutations",
//           [](const ARG& arg) {
//             vector<Mutation> mutations;
//             for (auto&& m : arg.mutations) {
//               mutations.push_back(*m);
//             }
//             return mutations;
//           },
//           "Return list of all Mutation objects")
//       .def(
//           "sites",
//           [](const ARG& arg) {
//             vector<Site> sites;
//             for (auto&& s : arg.sites) {
//               sites.push_back(*s);
//             }
//             return sites;
//           },
//           "Return list of all Site objects")
//       .def("num_nodes", &ARG::num_nodes)
//       .def("num_edges", &ARG::num_edges)
//       .def("num_mutations", &ARG::num_mutations)
//       .def(
//           "num_samples", [](const ARG& arg) { return arg.leaf_ids.size(); },
//           "Return number of samples")
//       .def_readonly("leaf_ids", &ARG::leaf_ids)         // unsorted
//       .def_readonly("sample_names", &ARG::sample_names) // unsorted
//       .def("set_offset", &ARG::set_offset, py::arg("offset"))
//       .def("set_chromosome", &ARG::set_chromosome, py::arg("chromosome"))
//       .def("set_sites", &ARG::set_sites, py::arg("positions"))
//       .def("add_sample", &ARG::add_sample, py::arg("sample_name") = "")
//       .def("is_leaf", &ARG::is_leaf, py::arg("node_id"))
//       .def("initialize_threading", &ARG::initialize_threading, py::arg("haps"))
//       .def("thread_sample",
//            py::overload_cast<vector<real>, vector<int>, vector<real>>(&ARG::thread_sample),
//            py::arg("section_starts"), py::arg("sample_ids"), py::arg("heights"))
//       .def("thread_sample",
//            py::overload_cast<vector<real>, vector<int>, vector<real>, vector<int>, vector<int>>(
//                &ARG::thread_sample),
//            py::arg("section_starts"), py::arg("sample_ids"), py::arg("heights"),
//            py::arg("thread_haps"), py::arg("target_haps"))
//       .def("populate_children_and_roots", &ARG::populate_children_and_roots)
//       .def("populate_mutations", &ARG::populate_mutations)
//       .def("lowest_mutated_edge", &ARG::lowest_mutated_edge, py::return_value_policy::reference)
//       .def("lowest_mutated_edge_by_site", &ARG::lowest_mutated_edge_by_site,
//            py::return_value_policy::reference)
//       .def("threading_constraints", &ARG::threading_constraints, py::arg("target_ID"), py::arg("g_thread"), py::arg("g_target"), py::arg("from"), py::arg("to"))
//       .def("root_starts", &ARG::root_starts)
//       .def("root_at", &ARG::root_at, py::return_value_policy::reference, py::arg("position"))
//       .def("mrca", &ARG::mrca, py::return_value_policy::reference, py::arg("ID1"), py::arg("ID2"),
//            py::arg("position"))
//       // overloaded case
//       // .def("mrca", (ARGNode* (ARG::*)(int, int, real)) &ARG::mrca,
//       //   py::arg("ID1"), py::arg("ID2"), py::arg("position"),
//       //   py::return_value_policy::reference)
//       // .def("mrca", (ARGNode* (ARG::*)(ARGNode&, ARGNode&, real)) &ARG::mrca,
//       //   py::arg("node1"), py::arg("node2"), py::arg("position"),
//       //   py::return_value_policy::reference)
//       .def("check_basic", &ARG::check_basic, "Basic checks after initializing an ARG",
//            py::arg("stringent") = true)
//       .def("check_roots", &ARG::check_roots, "Checks after populating roots")
//       .def("check_children", &ARG::check_children, "Checks after populating children",
//            py::arg("stringent") = true)
//       .def("__str__", [](const ARG& arg) {
//         std::ostringstream oss;
//         oss << arg;
//         return oss.str();
//       });

//   py::class_<DescendantList>(m, "DescendantList")
//       .def_static(
//           "set_threshold", &DescendantList::set_threshold, "Set threshold", py::arg("threshold"))
//       .def_static("print_threshold", &DescendantList::print_threshold, "Print threshold");

//   m.def("arg_to_newick", &arg_utils::arg_to_newick, py::arg("arg"), py::arg("verbose") = false,
//         "Return a Newick representation of an ARG, `verbose` includes branch lengths");
//   m.def("tmrca_mse", &arg_utils::tmrca_mse, py::arg("arg1"), py::arg("arg2"),
//         "Weighted TMRCA MSE between two ARGs");
//   m.def("kc_topology", &arg_utils::kc_topology, py::arg("arg1"), py::arg("arg2"),
//         "Weighted KC topology distance between two ARGs");
//   m.def("metrics_stab", &arg_utils::metrics_stab, py::arg("arg1"), py::arg("arg2"),
//         py::arg("num_stabs"),
//         "Weighted KC topology distance and TMRCA MSE between two ARGs, sampled based on stabbing "
//         "queries");
//   m.def("metrics_stab_efficient", &arg_utils::metrics_stab_efficient, py::arg("arg1"),
//         py::arg("arg2"), py::arg("num_stabs"), py::arg("random_kc_seed") = 0,
//         py::arg("merge_type") = 0, py::arg("merge_fraction") = 0, py::arg("use_r2") = false,
//         py::arg("use_log2") = false,
//         "Weighted KC topology distance and TMRCA MSE between two ARGs, sampled based on stabbing "
//         "queries");
//   m.def("bitset_overlap_full", &arg_utils::bitset_overlap_full, py::arg("arg1"), py::arg("arg2"),
//         py::arg("min_position") = -1, py::arg("max_position") = -1,
//         "Bitset overlap for branch precision and recall, over the full ARG or subregion");
//   m.def("bitset_overlap_stab", &arg_utils::bitset_overlap_stab, py::arg("arg1"), py::arg("arg2"),
//         py::arg("num_stabs"), py::arg("arg2_factor") = 1, py::arg("random_resolve_seed") = 0,
//         "Bitset overlap for branch and mutation precision and recall, sampled based on stabbing "
//         "queries");
//   m.def(
//       "stab_return_all_bitsets",
//       [](const ARG& arg, real position) {
//         vector<tuple<int, real, deque<bool>>> new_result;
//         for (auto& entry : arg_utils::stab_return_all_bitsets(arg, position)) {
//           new_result.emplace_back(
//               std::get<0>(entry), std::get<1>(entry), std::get<2>(entry).to_deque_bool());
//         }
//         return new_result;
//       },
//       py::arg("arg"), py::arg("position"),
//       "All bitsets at a position of the ARG, each tuple with values allele count / length / "
//       "bitset");
//   m.def("impute", &arg_utils::impute, py::arg("arg"), py::arg("position"), py::arg("genotypes"),
//         py::arg("old") = false, "ARG imputation");
//   m.def("mutation_match", &arg_utils::mutation_match, py::arg("arg"), py::arg("position"),
//         py::arg("genotypes"),
//         "Boolean for whether genotypes can arise from a single mutation at position");
//   m.def("mutation_best", &arg_utils::mutation_best, py::arg("arg"), py::arg("position"),
//         py::arg("genotypes"), py::arg("random_seed") = 0,
//         "Hamming distance for best mutation placement on ARG at a position");
//   m.def("distance_matrix", &arg_utils::distance_matrix, py::arg("arg"),
//         "Between-sample distance matrix in upper diagonal form");
//   m.def("distance_matrix_v2", &arg_utils::distance_matrix_v2, py::arg("arg"), py::arg("alpha") = 0,
//         "Between-sample distance matrix in upper diagonal form");
//   m.def("distance_matrix_time_bins", &arg_utils::distance_matrix_time_bins, py::arg("arg"),
//         py::arg("time_bins"),
//         "Between-sample distance matrices in upper diagonal form taking in time bins");
//   m.def("distance_matrix_maf_bins", &arg_utils::distance_matrix_maf_bins, py::arg("arg"),
//         py::arg("maf_bins"),
//         "Between-sample distance matrices in upper diagonal form taking in MAF bins");
//   m.def("association_diploid", &arg_utils::association_diploid, py::arg("arg"),
//         py::arg("raw_phenotypes"), py::arg("use_sample"), py::arg("file_root"),
//         py::arg("chromosome") = 1, py::arg("snp_prefix") = "",
//         py::arg("write_bitset_threshold") = -1, py::arg("calibration_factor") = 1,
//         py::arg("min_maf") = -1, py::arg("max_maf") = -1, py::arg("careful") = true,
//         py::arg("concise_pvalue") = true, py::arg("max_only") = false,
//         "ARG diploid branch association");
//   m.def("association_diploid_mutation", &arg_utils::association_diploid_mutation, py::arg("arg"),
//         py::arg("raw_phenotypes"), py::arg("use_sample"), py::arg("file_root"),
//         py::arg("chromosome") = 1, py::arg("snp_prefix") = "", py::arg("mu") = 0,
//         py::arg("random_seed") = 0, py::arg("write_bitset_threshold") = -1,
//         py::arg("calibration_factor") = 1, py::arg("min_maf") = -1, py::arg("max_maf") = -1,
//         py::arg("careful") = true, py::arg("concise_pvalue") = true, py::arg("max_only") = false,
//         "ARG diploid branch association");
//   m.def("association_diploid_mutation_multimu", &arg_utils::association_diploid_mutation_multimu,
//         py::arg("arg"), py::arg("raw_phenotypes"), py::arg("use_sample"), py::arg("file_root"),
//         py::arg("chromosome") = 1, py::arg("snp_prefix") = "", py::arg("mus") = vector<real>(0),
//         py::arg("random_seed") = 0, py::arg("write_bitset_threshold") = -1,
//         py::arg("calibration_factor") = 1, py::arg("min_maf") = -1, py::arg("max_maf") = -1,
//         py::arg("careful") = true, py::arg("concise_pvalue") = true, py::arg("max_only") = false,
//         "ARG diploid branch association");
//   m.def("association_max", &arg_utils::association_haploid_max, py::arg("arg"),
//         py::arg("raw_phenotypes"), py::arg("compress") = true, py::arg("epsilon") = 1e-8,
//         "ARG branch haploid association");
//   m.def("write_bitsets", &arg_utils::write_bitsets, py::arg("arg"), py::arg("file_root") = "",
//         py::arg("diploid") = false, py::arg("chromosome") = 1, py::arg("snp_prefix") = "",
//         py::arg("compress") = true, py::arg("count_only") = false,
//         "Write out ARG branches as bitsets");
//   m.def("write_bitsets_efficient", &arg_utils::write_bitsets_efficient, py::arg("arg"),
//         py::arg("file_root") = "", py::arg("diploid") = false, py::arg("chromosome") = 1,
//         py::arg("snp_prefix") = "", py::arg("min_mac") = 0, py::arg("max_mac") = 0,
//         py::arg("write_dosage") = false, py::arg("use_gz") = false, py::arg("count_only") = false,
//         "Write out ARG branches as bitsets");
//   m.def("visit_identical", &arg_utils::visit_identical, py::arg("arg"), py::arg("rel_tol") = 1e-6,
//         py::arg("abs_tol") = 0, py::arg("timing") = true, py::arg("verbose") = false,
//         "Check for identical ARG visit output between slow and fast versions");
//   m.def("time_efficient_visit", &arg_utils::time_efficient_visit, py::arg("arg"),
//         py::arg("timing") = false, "Time efficient visit routine");
//   m.def("pier_visit_arg", &arg_utils::visit_arg, py::arg("arg"), py::arg("verbose") = false,
//         "Pier's slow visit function which returns a map of bitstring to volume");
//   m.def("total_volume", &arg_utils::total_volume, py::arg("arg"), "Get ARG volume");
//   m.def("local_volume", &arg_utils::local_volume, py::arg("arg"), py::arg("min_position"),
//         py::arg("max_position"), "Get the local arg volume");
//   m.def("ARG_by_matrix_multiply_new_muts", &arg_utils::ARG_matrix_multiply_new_mut, py::arg("arg"),
//         py::arg("matrix"), py::arg("mu"), py::arg("random_seed") = 0, py::arg("min_time") = 0,
//         py::arg("max_time") = std::numeric_limits<double>::infinity(),
//         "Sample new mutations and multiply each by a matrix");
//   m.def("ARG_by_matrix_multiply_muts", &arg_utils::ARG_matrix_multiply_existing_mut, py::arg("arg"),
//         py::arg("matrix"), py::arg("standardize_mut") = false, py::arg("min_time") = 0,
//         py::arg("max_time") = std::numeric_limits<double>::infinity(),
//         "Multiply each existing ARG mutation by a sample-by-k matrix");
//   m.def("ARG_by_matrix_multiply_samples", &arg_utils::ARG_matrix_multiply_samples, py::arg("arg"),
//         py::arg("matrix"), py::arg("standardize_mut") = false, py::arg("min_time") = 0,
//         py::arg("max_time") = std::numeric_limits<double>::infinity(),
//         "Multiply each sample by a mutations-by-k matrix");
//   m.def("generate_mutations", &arg_utils::generate_mutations_and_return, py::arg("arg"),
//         py::arg("mu"), py::arg("random_seed") = 0, "Generate mutations and return a map");
//   m.def("write_mutations_to_haps", &arg_utils::write_mutations_to_haps, py::arg("arg"),
//         py::arg("file_root"), py::arg("min_maf") = 0, py::arg("max_maf") = 1,
//         py::arg("min_time") = 0, py::arg("max_time") = std::numeric_limits<double>::infinity(),
//         "Write mutations to disk in haps/samples format");
//   m.def("generate_m_mutations_and_keep", &arg_utils::generate_m_mutations_and_keep, py::arg("arg"),
//         py::arg("M"), py::arg("random_seed") = 0, py::arg("sorted") = true,
//         "Generate exactly M mutations and keep on the ARG");
//   m.def("generate_mutations_and_keep", &arg_utils::generate_mutations_and_keep, py::arg("arg"),
//         py::arg("mu"), py::arg("random_seed") = 0, py::arg("sorted") = true,
//         "Generate mutations with a given rate and keep on the ARG");
//   m.def("num_lineages", &arg_utils::num_lineages, py::arg("arg"), py::arg("position"),
//         py::arg("height"), "Count the number of lineages at a given position and height");
// }
