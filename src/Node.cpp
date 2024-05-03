#include "Node.hpp"

#include <iostream>


Node::Node(int _sample_ID, int _divergence, bool _genotype)
    : sample_ID(_sample_ID), divergence(_divergence), genotype(_genotype) {
}

void Node::insert_above(Node* node) {
  Node* old_above = above;
  old_above->below = node;
  above = node;
  node->above = old_above;
  node->below = this;
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
  os << "Node for sample " << node.sample_ID << " carrying allele " << node.genotype;
  return os;
}
