#include "State.hpp"

#include <memory>
#include <vector>



class DPBWT {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  

public:
  int num_sites;
  int num_samples;
  std::vector<int> physical_positions;
  std::vector<float> genetic_positions;
  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;
  float mutation_rate;

  // Constructors
  DPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate);

  // Insertion and extention functions
  void insert(std::vector<bool> genotype);
  Node* extend_node(Node* node, bool genotype);
  std::tuple<bool, int> extensible_by(State& s, const Node* t_next, const bool g);
  bool state_node_prefix_match(State state, int sample_ID);
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);
  // Algorithms
  std::vector<int> longest_prefix(std::vector<bool> genotype);
  std::vector<int> fastLS(std::vector<bool> genotype);
  double mutation_penalty();
  double recombination_penalty();

  // Debugging
  void print_sorting();
  void print_divergence();
};
