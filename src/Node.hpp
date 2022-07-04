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

  // "Next below to the right" for 0 and 1
  std::vector<Node*> w = std::vector<Node*>(2);

  // Constructors
  Node(int sample_ID, int site, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);

  // Output
  friend ostream& operator<<(ostream& os, const Node& node);
};