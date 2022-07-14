#include "Node.hpp"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

// Node::Node() {}

// Node::Node(int _sample_ID, int _site, bool _genotype)
//   : sample_ID(_sample_ID), site(_site), genotype(_genotype),
Node::Node(int _sample_ID, bool _genotype)
  : sample_ID(_sample_ID), genotype(_genotype),
  above(nullptr),
  below(nullptr),
  w({nullptr, nullptr}) {
  // divergence(_site) {
}

void Node::insert_above(Node* node) {
	Node* old_above = above;
  old_above->below = node;
	above = node;
	node->above = old_above;
	node->below = this;
}

ostream& operator<<(ostream& os, const Node& node) {
  os << "Node for sample " << node.sample_ID << " carrying allele " << node.genotype;
  return os;
}
