#include "State.hpp"
#include "HMM.hpp"

#include <memory>
#include <vector>
#include <unordered_map>


class Threads {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;  // Ha, ha

  Node* extend_node(Node* node, bool genotype, const int i);
  bool extensible_by(State& s, const Node* t_next, const bool g, const int i);
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);
  std::vector<bool> fetch_het_hom_sites(const int id1, const int id2, const int start, const int end);
  std::vector<int> het_sites_from_thread(const int focal_ID, std::vector<int> bp_starts, std::vector<int> target_IDs);

public:
  int n_prune;
  int num_sites;
  int num_samples;
  double mutation_rate;
  // This determines the segment dating formula
  bool sparse_sites;
  // Whether to internally break up big segments
  bool use_hmm;
  std::unordered_map<int, int> ID_map;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  std::vector<double> bp_boundaries;
  std::vector<double> cm_boundaries;
  Demography demography;
  HMM* hmm;

  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;

  // Constructors and utils
  Threads(std::vector<double> _physical_positions, std::vector<double> _genetic_positions, double _mutation_rate, double ne, bool _sparse_sites, int _n_prune, bool _use_hmm) :
    Threads(_physical_positions, _genetic_positions, _mutation_rate, std::vector<double>{ne}, std::vector<double>{0.0}, _sparse_sites, _n_prune, _use_hmm) {};
  Threads(std::vector<double> _physical_positions, std::vector<double> _genetic_positions, double _mutation_rate, std::vector<double> ne, std::vector<double> ne_times, bool _sparse_sites, int _n_prune, bool _use_hmm);
  std::tuple<std::vector<double>, std::vector<double>> site_sizes(std::vector<double> positions);

  // Insertion/deletion
  // Insert and assign generic ID
  void insert(const std::vector<bool>& genotype);
  // Insert and assign specific ID
  void insert(const int ID, const std::vector<bool>& genotype);
  // Remove sample from panel
  void delete_ID(int ID);

  // HMM
  void delete_hmm();

  // Algorithms
  std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, std::vector<int>> thread(const std::vector<bool>& genotype);
  std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, std::vector<int>> thread(const int new_sample_ID, const std::vector<bool>& genotype);

  // std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, std::vector<int>, std::vector<bool>> thread_with_mutations(const std::vector<bool>& genotype);
  // std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, std::vector<int>, std::vector<bool>> thread_with_mutations(const int new_sample_ID, const std::vector<bool>& genotype);
  // std::vector<std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>> thread_from_file(std::string file_path, int n_cycle);
  std::vector<std::tuple<int, int>> fastLS(const std::vector<bool>& genotype);
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties();
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties_correct();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties_correct();

  // TMRCA estimation
  double date_segment(const int id1, const int id2, const int start, const int end);

  // Debugging
  void print_sorting();
};
