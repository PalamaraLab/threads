// This file is part of the Threads software suite.
// Copyright (C) 2024-2025 Threads Developers.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
