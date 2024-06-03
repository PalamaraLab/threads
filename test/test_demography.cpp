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

#include "Demography.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

TEST_CASE("Demography constructor exceptions") {
  SECTION("Throws if times and sizes vector lengths differ") {
    CHECK_THROWS_WITH(Demography({0, 0, 0}, {0}), Catch::Matchers::ContainsSubstring("Demography times and sizes must have equal length"));
    CHECK_THROWS_WITH(Demography({0}, {0, 0, 0}), Catch::Matchers::ContainsSubstring("Demography times and sizes must have equal length"));
  }

  SECTION("Throws if times not non-negative") {
    CHECK_THROWS_WITH(Demography({1}, {-1}), Catch::Matchers::ContainsSubstring("Demography expects non-negative times"));
  }

  SECTION("Throws if sizes not positive") {
    CHECK_THROWS_WITH(Demography({0}, {0}), Catch::Matchers::ContainsSubstring("strictly positive sizes"));
    CHECK_THROWS_WITH(Demography({-1}, {0}), Catch::Matchers::ContainsSubstring("strictly positive sizes"));
  }

  SECTION("Throws if times not increasing") {
    CHECK_THROWS_WITH(Demography({1, 1, 1}, {3, 2, 1}), Catch::Matchers::ContainsSubstring("Demography times must be strictly increasing. Found 2 after 3 at index 1"));
  }

  SECTION("Throws if times does not start at 0") {
    CHECK_THROWS_WITH(Demography({1}, {5}), Catch::Matchers::ContainsSubstring("Demography must start at time 0.0"));
  }
}

TEST_CASE("Demography std_to_gen exceptions") {
  SECTION("std_to_gen value must be non-negative") {
    Demography d{{1}, {0}};
    CHECK_THROWS_WITH(d.std_to_gen(-1), Catch::Matchers::ContainsSubstring("Demography can only convert times greater than the first entry"));
  }
}
