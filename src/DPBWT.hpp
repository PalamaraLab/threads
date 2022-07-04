#include "Node.hpp"

#include <memory>
#include <vector>

class State {
public:
  Node* below;
  double score;
  State* parent;
  int start;
  bool genotype;
  bool extended_by_div;
  int best_above;

  State(Node* _below, double _score, State* _parent, int _start, bool _genotype, bool _extended_by_div, int _best_above);

  friend ostream& operator<<(ostream& os, State& state);
};

class DPBWT {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  // The dynamic reference panel
  // To keep stack-states alive during fastLS inference
  std::vector<std::unique_ptr<State>> states;

public:
  int num_sites;
  int num_samples;
  std::vector<int> physical_positions;
  std::vector<float> genetic_positions;
  std::vector<std::vector<std::unique_ptr<Node>>> panel;
  float mutation_rate;

  // Constructors
  DPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate);

  // Insertion and extention functions
  void insert(std::vector<bool> genotype);
  Node* extend_node(Node* node, bool genotype);
  std::tuple<bool, bool, int> extensible_by(State s, Node* t_next, bool g);
  bool state_node_prefix_match(State state, int sample_ID);

  // Algorithms
  std::vector<int> longest_prefix(std::vector<bool> genotype);
  std::vector<int> fastLS(std::vector<bool> genotype);
  double mutation_penalty();
  double recombination_penalty();

  // Debugging
  void print_sorting();
  void print_divergence();
};
