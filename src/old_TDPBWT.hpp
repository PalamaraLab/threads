#include "Node.hpp"

// #include <pair>
#include <memory>
#include <vector>
// using std::vector;

class TDPBWT {

// private:
  // void _make_mismatch_probs();
  // bool feasible();
  // This is the back of each (horizontal) deque. Dummy nodes.
  // This is the top of each (vertical) deque. Actual nodes.
  std::vector<std::unique_ptr<Node>> all_nodes;
  std::vector<Node*> backs;
  std::vector<Node*> tops;
  std::vector<Node*> bottoms;
  Node virtual_top;
  Node virtual_bottom;

public:
  int num_sites;
  int num_samples;
  std::vector<int> physical_positions;
  std::vector<float> genetic_positions;
  float mutation_rate;

  TDPBWT(std::vector<int> _physical_positions, std::vector<float> _genetic_positions, float _mutation_rate);

  void insert(std::vector<bool> genotype);
  // std::vector<int> best_path_haploid(std::vector<bool> genotype);
  // std::vector<int> best_path_diploid(std::vector<bool> genotype);
  std::pair<int, int> longest_suffix(std::vector<bool> genotype);
  void print_sorting();
  void print_U0();
  void print_U1();
};
