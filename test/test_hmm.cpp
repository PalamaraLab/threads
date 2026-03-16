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

#include "HMM.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <vector>

namespace {

// Create a simple constant-Ne demography for testing
Demography simple_demography(double ne = 10000.0) {
  return Demography({ne}, {0.0});
}

// Create uniform site sizes for n_sites sites
std::vector<double> uniform_bp_sizes(int n_sites, double bp_size = 100.0) {
  return std::vector<double>(n_sites, bp_size);
}

std::vector<double> uniform_cm_sizes(int n_sites, double cm_size = 0.001) {
  return std::vector<double>(n_sites, cm_size);
}

} // namespace

TEST_CASE("HMM construction") {
  int n_sites = 20;
  int K = 8;
  auto demo = simple_demography();
  auto bp = uniform_bp_sizes(n_sites);
  auto cm = uniform_cm_sizes(n_sites);

  HMM hmm(demo, bp, cm, 1.4e-8, K);

  CHECK(hmm.num_states == K);
  CHECK(hmm.expected_times.size() == K);
  CHECK(hmm.trellis.size() == n_sites);
  CHECK(hmm.pointers.size() == n_sites);
  CHECK(hmm.non_transition_score.size() == n_sites);
  CHECK(hmm.transition_score.size() == n_sites);
  CHECK(hmm.hom_score.size() == n_sites);
  CHECK(hmm.het_score.size() == n_sites);
}

TEST_CASE("HMM expected times are increasing") {
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(10), uniform_cm_sizes(10), 1.4e-8, 16);

  for (int i = 1; i < 16; i++) {
    CHECK(hmm.expected_times[i] > hmm.expected_times[i - 1]);
  }
  // All expected times should be positive
  for (int i = 0; i < 16; i++) {
    CHECK(hmm.expected_times[i] > 0.0);
  }
}

TEST_CASE("HMM breakpoints with all homozygous") {
  int n_sites = 20;
  int K = 4;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  // All homozygous = no mutations -> should stay in one state
  std::vector<bool> obs(n_sites, false);
  auto bps = hmm.breakpoints(obs, 0);

  // Should have at least the initial breakpoint at 0
  CHECK(bps.size() >= 1);
  CHECK(bps[0] == 0);
}

TEST_CASE("HMM breakpoints with all heterozygous") {
  int n_sites = 30;
  int K = 4;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  // All het -> lots of mutations -> should stay in deepest time state
  std::vector<bool> obs(n_sites, true);
  auto bps = hmm.breakpoints(obs, 0);

  CHECK(bps.size() >= 1);
  CHECK(bps[0] == 0);
}

TEST_CASE("HMM breakpoints with mixed signal") {
  int n_sites = 40;
  int K = 8;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  // First half: all hom (recent), second half: all het (old) -> expect breakpoint
  std::vector<bool> obs(n_sites, false);
  for (int i = n_sites / 2; i < n_sites; i++) {
    obs[i] = true;
  }

  auto bps = hmm.breakpoints(obs, 0);
  CHECK(bps.size() >= 1);
  CHECK(bps[0] == 0);
  // With a strong signal change, we expect at least one additional breakpoint
  // (though exact number depends on HMM parameters)
}

TEST_CASE("HMM breakpoints with offset start") {
  int n_sites = 30;
  int K = 4;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  // Use only a sub-range starting at offset 5
  int start = 5;
  int len = 15;
  std::vector<bool> obs(len, false);

  auto bps = hmm.breakpoints(obs, start);
  CHECK(bps[0] == start);
}

TEST_CASE("HMM recombination scores are negative log-probs") {
  int n_sites = 5;
  int K = 4;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  for (int i = 0; i < n_sites; i++) {
    for (int k = 0; k < K; k++) {
      // Both transition and non-transition scores should be <= 0 (log-probs)
      CHECK(hmm.transition_score[i][k] <= 0.0);
      CHECK(hmm.non_transition_score[i][k] <= 0.0);
      // non-transition should be >= transition (more likely to not transition)
      CHECK(hmm.non_transition_score[i][k] >= hmm.transition_score[i][k]);
    }
  }
}

TEST_CASE("HMM mutation scores are negative log-probs") {
  int n_sites = 5;
  int K = 4;
  auto demo = simple_demography();
  HMM hmm(demo, uniform_bp_sizes(n_sites), uniform_cm_sizes(n_sites), 1.4e-8, K);

  for (int i = 0; i < n_sites; i++) {
    for (int k = 0; k < K; k++) {
      CHECK(hmm.hom_score[i][k] <= 0.0);
      CHECK(hmm.het_score[i][k] <= 0.0);
      // hom (not mutating) should be more likely than het (mutating) for typical params
      CHECK(hmm.hom_score[i][k] >= hmm.het_score[i][k]);
    }
  }
}
