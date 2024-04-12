#ifndef THREADS_INFER_NODE_HPP
#define THREADS_INFER_NODE_HPP

#include <array>
#include <iostream>

using std::ostream;

class Node {
public:
  // Node data
  int sample_ID;
  int divergence;
  bool genotype;

  // Linked list pointers
  Node* above = nullptr;
  Node* below = nullptr;

  // "Next below to the right" for 0 and 1
  std::array<Node*, 2> w = {nullptr, nullptr};

  // Constructors
  Node(int sample_ID, int divergence, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);

  // Output
  friend ostream& operator<<(ostream& os, const Node& node);
};

#endif // THREADS_INFER_NODE_HPP