#ifndef THREADS_ARG_NODE_HPP
#define THREADS_ARG_NODE_HPP

#include <array>
#include <iostream>

class Node {
public:
  // Constructors
  Node(int sample_ID, int divergence, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);

  // Output
  friend std::ostream& operator<<(std::ostream& os, const Node& node);

public:
  // Node data
  int sample_ID = 0;
  int divergence = 0;
  bool genotype = 0;

  // Linked list pointers
  Node* above = nullptr;
  Node* below = nullptr;

  // "Next below to the right" for 0 and 1
  std::array<Node*, 2> w = {nullptr, nullptr};
};

#endif // THREADS_ARG_NODE_HPP