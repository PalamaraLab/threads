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

#include <catch2/catch_test_macros.hpp>
#include <sstream>

TEST_CASE("Node construction") {
  Node n(42, 10, true);
  CHECK(n.sample_ID == 42);
  CHECK(n.divergence == 10);
  CHECK(n.genotype == true);
  CHECK(n.above == nullptr);
  CHECK(n.below == nullptr);
  CHECK(n.w[0] == nullptr);
  CHECK(n.w[1] == nullptr);
}

TEST_CASE("Node insert_above") {
  // Set up a two-node chain: bottom <-> top
  Node bottom(0, 0, false);
  Node top(1, 0, true);
  bottom.above = &top;
  top.below = &bottom;

  // Insert middle between bottom and top
  Node middle(2, 5, false);
  bottom.insert_above(&middle);

  // Verify chain is now bottom <-> middle <-> top
  CHECK(bottom.above == &middle);
  CHECK(middle.below == &bottom);
  CHECK(middle.above == &top);
  CHECK(top.below == &middle);
}

TEST_CASE("Node insert_above multiple") {
  // Build chain of 4 nodes by inserting above bottom
  Node bottom(0, 0, false);
  Node top(1, 0, true);
  bottom.above = &top;
  top.below = &bottom;

  Node n1(2, 1, true);
  Node n2(3, 2, false);

  bottom.insert_above(&n1);
  n1.insert_above(&n2);

  // Chain: bottom <-> n1 <-> n2 <-> top
  CHECK(bottom.above == &n1);
  CHECK(n1.above == &n2);
  CHECK(n2.above == &top);
  CHECK(top.below == &n2);
}

TEST_CASE("Node stream output") {
  Node n(7, 3, true);
  std::ostringstream oss;
  oss << n;
  CHECK(oss.str() == "Node for sample 7 carrying allele 1");

  Node n2(0, 0, false);
  std::ostringstream oss2;
  oss2 << n2;
  CHECK(oss2.str() == "Node for sample 0 carrying allele 0");
}
