#include "State.hpp"

#include <memory>
#include <vector>



class DPBWT {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;

public:
  int num_sites;
  int num_samples;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;
  double mutation_rate;

  // Constructors and utils
  DPBWT(std::vector<double> _physical_positions, std::vector<double> _genetic_positions, double _mutation_rate);
  std::vector<double> site_sizes(std::vector<double> positions);

  // Insertion and extention functions
  void insert(std::vector<bool> genotype);
  Node* extend_node(Node* node, bool genotype);
  std::tuple<bool, int> extensible_by(State& s, const Node* t_next, const bool g);
  bool state_node_prefix_match(State state, int sample_ID);
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);

  // Algorithms
  std::vector<int> longest_prefix(std::vector<bool> genotype);
  std::vector<int> fastLS(std::vector<bool> genotype);
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties();

  // Debugging
  void print_sorting();
  void print_divergence();
};
