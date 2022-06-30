#include "Gnode.hpp"

// #include <pair>
#include <memory>
#include <vector>
// using std::vector;

class DPBWT {

private:
  // void _make_mismatch_probs();
  // bool feasible();
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  Node virtual_top;
  Node virtual_bottom;
  std::vector<std::vector<std::unique_ptr<Node>>> panel;

public:
  int num_sites;
  int num_samples;
  std::vector<int> physical_positions;
  std::vector<float> genetic_positions;
  float mutation_rate;

  DPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate);

  void insert(std::vector<bool> genotype);
  // std::pair<int, int> longest_prefix(std::vector<bool> genotype);
  void print_sorting();
  // void print_U0();
  // void print_U1();
};
