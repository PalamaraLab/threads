#include "Gnode.hpp"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

// Node::Node() {}

Node::Node(int _sample_ID, int _site, bool _genotype)
  : sample_ID(_sample_ID), site(_site), genotype(_genotype),
  above(nullptr),
  below(nullptr),
  left(nullptr),
  right(nullptr),
  w({nullptr, nullptr}) {
}

void Node::insert_above(Node* node) {
	Node* old_above = above;
  old_above->below = node;
	above = node;
	node->above = old_above;
	node->below = this;
  cout << "insert done\n";
}

bool Node::equivalent_to(Node* other) {
  return true;
    // return next_1 == other->next_1 && next_0 == other->next_0 && prev_1 == other->prev_1 && prev_0 == other->prev_0;
}

ostream& operator<<(ostream& os, const Node& node) {
  os << "Node for sample " << node.sample_ID << " and site " << node.site;
  return os;
}