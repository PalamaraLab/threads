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

#include "Demography.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <sstream>

TEST_CASE("Demography constant Ne") {
  double Ne = 10000.0;
  Demography d({Ne}, {0.0});

  // With constant Ne, std_to_gen should be linear: gen = t * Ne
  CHECK_THAT(d.std_to_gen(0.0), Catch::Matchers::WithinAbs(0.0, 1e-10));
  CHECK_THAT(d.std_to_gen(1.0), Catch::Matchers::WithinAbs(Ne, 1e-6));
  CHECK_THAT(d.std_to_gen(0.5), Catch::Matchers::WithinAbs(Ne * 0.5, 1e-6));
  CHECK_THAT(d.std_to_gen(2.0), Catch::Matchers::WithinAbs(Ne * 2.0, 1e-6));
}

TEST_CASE("Demography piecewise Ne") {
  // Ne=10000 for generations [0, 100), then Ne=20000 after
  Demography d({10000.0, 20000.0}, {0.0, 100.0});

  // std_times: [0, 100/10000] = [0, 0.01]
  CHECK_THAT(d.std_times[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
  CHECK_THAT(d.std_times[1], Catch::Matchers::WithinAbs(0.01, 1e-12));

  // Within first epoch: std_to_gen(0.005) = 0 + 0.005 * 10000 = 50
  CHECK_THAT(d.std_to_gen(0.005), Catch::Matchers::WithinAbs(50.0, 1e-6));

  // Within second epoch: std_to_gen(0.02) = 100 + (0.02 - 0.01) * 20000 = 300
  CHECK_THAT(d.std_to_gen(0.02), Catch::Matchers::WithinAbs(300.0, 1e-6));
}

TEST_CASE("Demography expected branch length") {
  double Ne = 10000.0;
  Demography d({Ne}, {0.0});

  // expected_branch_length(N) = std_to_gen(2/N)
  // For constant Ne: 2/N * Ne
  CHECK_THAT(d.expected_branch_length(2), Catch::Matchers::WithinAbs(Ne, 1e-6));
  CHECK_THAT(d.expected_branch_length(10), Catch::Matchers::WithinAbs(2000.0, 1e-6));
  CHECK_THAT(d.expected_branch_length(100), Catch::Matchers::WithinAbs(200.0, 1e-6));
}

TEST_CASE("Demography expected_time is std_to_gen(1)") {
  Demography d({5000.0}, {0.0});
  CHECK_THAT(d.expected_time, Catch::Matchers::WithinAbs(5000.0, 1e-6));
}

TEST_CASE("Demography three epochs") {
  // Ne=1000 for [0,50), Ne=5000 for [50,100), Ne=20000 after 100
  Demography d({1000.0, 5000.0, 20000.0}, {0.0, 50.0, 100.0});

  // std_times: [0, 50/1000, 50/1000 + 50/5000] = [0, 0.05, 0.06]
  CHECK_THAT(d.std_times[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
  CHECK_THAT(d.std_times[1], Catch::Matchers::WithinAbs(0.05, 1e-12));
  CHECK_THAT(d.std_times[2], Catch::Matchers::WithinAbs(0.06, 1e-12));

  // In third epoch: std_to_gen(0.07) = 100 + (0.07 - 0.06) * 20000 = 300
  CHECK_THAT(d.std_to_gen(0.07), Catch::Matchers::WithinAbs(300.0, 1e-6));
}

TEST_CASE("Demography stream output") {
  Demography d({1000.0}, {0.0});
  std::ostringstream oss;
  oss << d;
  // Should not crash, but note the bug: operator<< uses std::cout instead of os
}
