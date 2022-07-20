#include "State.hpp"

#include <memory>
#include <vector>
#include <unordered_map>


class DPBWT {

private:
  // To access the linked lists in each column
  std::vector<std::unique_ptr<Node>> tops;
  std::vector<std::unique_ptr<Node>> bottoms;
  std::unordered_map<int, int> ID_map;

  Node* extend_node(Node* node, bool genotype, const int i);
  bool extensible_by(State& s, const Node* t_next, const bool g, const int i);
  bool genotype_interval_match(const int id1, const int id2, const int start, const int end);

public:
  double ne;
  int num_sites;
  int num_samples;
  std::vector<double> physical_positions;
  std::vector<double> genetic_positions;
  std::vector<double> bp_sizes;
  std::vector<double> cm_sizes;
  // The dynamic reference panel
  std::vector<std::vector<std::unique_ptr<Node>>> panel;
  double mutation_rate;

  // Constructors and utils
  DPBWT(std::vector<double> _physical_positions, std::vector<double> _genetic_positions, double _mutation_rate, double _ne);
  std::vector<double> site_sizes(std::vector<double> positions);

  // Insertion/deletion
  void insert(std::vector<bool> genotype);
  void insert(std::vector<bool> genotype, int ID);
  void delete_ID(int ID);

  // Algorithms
  std::vector<std::tuple<int, int, double, double>> thread(std::vector<bool> genotype);
  std::vector<int> longest_prefix(std::vector<bool> genotype);
  std::vector<std::tuple<int, int>> fastLS(std::vector<bool> genotype);
  std::tuple<std::vector<double>, std::vector<double>> mutation_penalties();
  std::tuple<std::vector<double>, std::vector<double>> recombination_penalties();

  // Aging
  double age_segment_ML(const int id1, const int id2, const int start, const int end);
  // (need to fix this)
  double age_segment_bayes(const int id1, const int id2, const int start, const int end);

  // Debugging
  void print_sorting();
};
