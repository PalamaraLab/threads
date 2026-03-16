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

#include "ThreadingInstructions.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <limits>
#include <vector>

TEST_CASE("ThreadingInstruction construction") {
  std::vector<int> starts = {100, 500, 800};
  std::vector<double> tmrcas = {50.0, 200.0, 100.0};
  std::vector<int> targets = {0, 3, 1};
  std::vector<int> mismatches = {2, 7, 15};

  ThreadingInstruction ti(starts, tmrcas, targets, mismatches);
  CHECK(ti.num_segments == 3);
  CHECK(ti.num_mismatches == 3);
  CHECK(ti.starts == starts);
  CHECK(ti.tmrcas == tmrcas);
  CHECK(ti.targets == targets);
  CHECK(ti.mismatches == mismatches);
}

TEST_CASE("ThreadingInstruction mismatching lengths throw") {
  CHECK_THROWS(ThreadingInstruction({0}, {1.0, 2.0}, {0}, {}));
  CHECK_THROWS(ThreadingInstruction({0, 1}, {1.0}, {0, 1}, {}));
}

TEST_CASE("ThreadingInstructions construction from components") {
  std::vector<std::vector<int>> starts = {{100, 500}, {100, 300, 700}};
  std::vector<std::vector<double>> tmrcas = {{10.0, 20.0}, {5.0, 15.0, 25.0}};
  std::vector<std::vector<int>> targets = {{0, 1}, {0, 1, 0}};
  std::vector<std::vector<int>> mismatches = {{1}, {0, 3}};
  std::vector<int> positions = {100, 200, 300, 400, 500, 600, 700, 800};

  ThreadingInstructions ti(starts, tmrcas, targets, mismatches, positions, 100, 800);
  CHECK(ti.num_samples == 2);
  CHECK(ti.num_sites == 8);
  CHECK(ti.start == 100);
  CHECK(ti.end == 800);
}

TEST_CASE("ThreadingInstructions all_starts/tmrcas/targets/mismatches") {
  std::vector<std::vector<int>> starts = {{0, 5}};
  std::vector<std::vector<double>> tmrcas = {{10.0, 20.0}};
  std::vector<std::vector<int>> targets = {{3, 7}};
  std::vector<std::vector<int>> mismatches = {{2}};
  std::vector<int> positions = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

  ThreadingInstructions ti(starts, tmrcas, targets, mismatches, positions, 100, 1000);

  auto all_s = ti.all_starts();
  CHECK(all_s.size() == 1);
  CHECK(all_s[0] == std::vector<int>{0, 5});

  auto all_t = ti.all_tmrcas();
  CHECK(all_t.size() == 1);
  CHECK(all_t[0][0] == 10.0);

  auto all_tg = ti.all_targets();
  CHECK(all_tg[0] == std::vector<int>{3, 7});

  auto all_m = ti.all_mismatches();
  CHECK(all_m[0] == std::vector<int>{2});
}

TEST_CASE("ThreadingInstructionIterator basic iteration") {
  std::vector<int> starts = {100, 500};
  std::vector<double> tmrcas = {10.0, 20.0};
  std::vector<int> targets = {3, 7};
  std::vector<int> mismatches = {2}; // mismatch at position index 2

  ThreadingInstruction ti(starts, tmrcas, targets, mismatches);
  std::vector<int> positions = {100, 200, 300, 400, 500, 600};

  ThreadingInstructionIterator iter(ti, positions);
  CHECK(iter.current_target == 3);
  CHECK(iter.current_tmrca == 10.0);

  // Advance past second segment start
  iter.increment_site(500);
  CHECK(iter.current_target == 7);
  CHECK(iter.current_tmrca == 20.0);
}

TEST_CASE("ThreadingInstructionIterator mismatch tracking") {
  std::vector<int> starts = {100};
  std::vector<double> tmrcas = {10.0};
  std::vector<int> targets = {0};
  std::vector<int> mismatches = {2}; // mismatch at site index 2

  ThreadingInstruction ti(starts, tmrcas, targets, mismatches);
  std::vector<int> positions = {100, 200, 300, 400, 500};

  ThreadingInstructionIterator iter(ti, positions);

  iter.increment_site(100);
  CHECK(iter.is_mismatch == false);

  iter.increment_site(200);
  CHECK(iter.is_mismatch == false);

  // Position 300 = positions[2] = the mismatch site
  iter.increment_site(300);
  CHECK(iter.is_mismatch == true);

  iter.increment_site(400);
  CHECK(iter.is_mismatch == false);
}
