#include <iostream>
#include <vector>
using std::ostream;

class Node {
public:
  // Node data
  int sample_ID;
  int site;
  bool genotype;
  int divergence;

  // Node pointers
  Node* above;
  Node* below;
  Node* left;
  Node* right;
  // "Next below to the right" for 0 and 1
  std::vector<Node*> w = std::vector<Node*>(2);

  // Constructors
  // Node();
  Node(Node&&) = default;
  Node(int sample_ID, int site, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);

  // Comparators
  // bool equivalent_to(Node* other);
  friend bool operator==(const Node& node1, const Node& node2) {
    return (node1.sample_ID == node2.sample_ID && node1.site == node2.site) && node1.genotype == node2.genotype;
  }

  // Output
  friend ostream& operator<<(ostream& os, const Node& node);
};