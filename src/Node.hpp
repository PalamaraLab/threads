// This file is part of the Threads software suite.
// Copyright (C) 2024 Threads Developers.
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

#ifndef THREADS_ARG_NODE_HPP
#define THREADS_ARG_NODE_HPP

#include <array>
#include <iostream>

class Node {
public:
  // Constructors
  Node() = default;
  Node(int sample_ID, int divergence, bool genotype);

  // Node movers-arounders
  void insert_above(Node* node);
  void assign(int sample_ID, int divergence, bool genotype);

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