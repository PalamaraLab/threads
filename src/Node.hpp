#include <vector>
#include <iostream>
using std::ostream;

class Node {
public:
  // Node data
  int sample_ID;
  bool genotype;

  // Linked list pointers
  Node* above;
  Node* below;

  // "Next below to the right" for 0 and 1
  Node* w[2];

  // Constructors
  Node(int sample_ID, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);

  // Output
  friend ostream& operator<<(ostream& os, const Node& node);
};
