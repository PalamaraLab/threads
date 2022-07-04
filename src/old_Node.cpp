#include "Node.hpp"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

Node::Node(int _sample_ID, int _site, bool _genotype)
  : sample_ID(_sample_ID), site(_site), genotype(_genotype),
  above(nullptr),
  below(nullptr),
  left(nullptr),
  right(nullptr),
  next_0(nullptr),
  next_1(nullptr),
  prev_0(nullptr),
  prev_1(nullptr) {
}

Node* Node::first_uncle(bool g) {
  if (above == nullptr) {
    return nullptr;
  } else if (above->left->genotype == g) {
    return above->left;
  } else {
    return above->first_uncle(g);
  }
}

Node* Node::first_nephew(bool g) {
  if (below == nullptr) {
    return nullptr;
  } else if (below->left->genotype == g) {
    return below->left;
  } else {
    return below->first_nephew(g);
  }
}

void Node::insert_right(Node* node) {
	Node* old_right = right;
  if (old_right != nullptr) {
    old_right->left = node;
  }
	right = node;
	node->right = old_right;
	node->left = this;
}

void Node::insert_left(Node* node) {
	Node* old_left = left;
	left = node;
	node->left = old_left;
	node->right = this;
}

void Node::insert_below(Node* node) {
	Node* old_below = below;
  if (old_below != nullptr) {
    old_below->above = node;
  }
	below = node;
	node->below = old_below;
	node->above = this;
}

void Node::insert_above(Node* node) {
	Node* old_above = above;
  if (old_above != nullptr) {
    old_above->below = node;
  }
	above = node;
	node->above = old_above;
	node->below = this;
}

void Node::update_next_0_recursive(Node* old_node, Node* new_node) {
	if (*next_0 == *old_node) {
		next_0 = new_node;
		if (above != nullptr) {
			above->update_next_0_recursive(old_node, new_node);
		}
	}
}

void Node::update_prev_0_recursive(Node* old_node, Node* new_node) {
  if (*prev_0 == *old_node) {
    prev_0 = new_node;
    if (below != nullptr) {
      below->update_prev_0_recursive(old_node, new_node);
    }
  }
}

void Node::update_next_1_recursive(Node* old_node, Node* new_node) {
  if (*next_1 == *old_node) {
    next_1 = new_node;
    if (above != nullptr) {
      above->update_next_1_recursive(old_node, new_node);
    }
  }
}

void Node::update_prev_1_recursive(Node* old_node, Node* new_node) {
  if (*prev_1 == *old_node) {
    prev_1 = new_node;
    if (below != nullptr) {
      below->update_prev_1_recursive(old_node, new_node);
    }
  }
}

bool Node::equivalent_to(Node* other) {
    return next_1 == other->next_1 && next_0 == other->next_0 && prev_1 == other->prev_1 && prev_0 == other->prev_0;
}

ostream& operator<<(ostream& os, const Node& node) {
  os << "Node for sample " << node.sample_ID << " and site " << node.site;
  return os;
}