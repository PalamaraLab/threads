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

#include "ViterbiLowMem.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <vector>

// === ViterbiPath tests ===

TEST_CASE("ViterbiPath construction and basic operations") {
  ViterbiPath path(5);
  CHECK(path.target_id == 5);
  CHECK(path.size() == 0);

  path.append(0, 3);
  path.append(10, 7);
  CHECK(path.size() == 2);
  CHECK(path.segment_starts[0] == 0);
  CHECK(path.segment_starts[1] == 10);
  CHECK(path.sample_ids[0] == 3);
  CHECK(path.sample_ids[1] == 7);
}

TEST_CASE("ViterbiPath reverse") {
  ViterbiPath path(0);
  path.append(0, 1);
  path.append(5, 2);
  path.append(10, 3);

  path.reverse();
  CHECK(path.segment_starts[0] == 10);
  CHECK(path.segment_starts[1] == 5);
  CHECK(path.segment_starts[2] == 0);
  CHECK(path.sample_ids[0] == 3);
  CHECK(path.sample_ids[1] == 2);
  CHECK(path.sample_ids[2] == 1);
}

TEST_CASE("ViterbiPath append with height and het_sites") {
  ViterbiPath path(0);
  std::vector<int> hets1 = {2, 3};
  path.append(0, 1, 100.0, hets1);
  std::vector<int> hets2 = {7};
  path.append(5, 2, 200.0, hets2);

  CHECK(path.size() == 2);
  CHECK(path.heights[0] == 100.0);
  CHECK(path.heights[1] == 200.0);
  CHECK(path.het_sites.size() == 3);
  CHECK(path.het_sites[0] == 2);
  CHECK(path.het_sites[1] == 3);
  CHECK(path.het_sites[2] == 7);
}

TEST_CASE("ViterbiPath append validates ordering") {
  ViterbiPath path(0);
  std::vector<int> hets = {};
  path.append(10, 1, 100.0, hets);

  // Appending segment_start <= previous should throw
  CHECK_THROWS(path.append(10, 2, 200.0, hets));
  CHECK_THROWS(path.append(5, 2, 200.0, hets));
}

TEST_CASE("ViterbiPath map_positions") {
  ViterbiPath path(0);
  path.append(0, 1);
  path.append(3, 2);
  path.append(7, 3);

  std::vector<int> positions = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  path.map_positions(positions);

  CHECK(path.bp_starts.size() == 3);
  CHECK(path.bp_starts[0] == 100);
  CHECK(path.bp_starts[1] == 400);
  CHECK(path.bp_starts[2] == 800);
}

TEST_CASE("ViterbiPath dump_data_in_range full") {
  ViterbiPath path(0);
  std::vector<int> hets1 = {};
  path.append(0, 1, 10.0, hets1);
  path.append(5, 2, 20.0, hets1);
  path.append(10, 3, 30.0, hets1);

  std::vector<int> positions = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100};
  path.map_positions(positions);

  auto [starts, ids, heights] = path.dump_data_in_range(-1, -1);
  CHECK(starts.size() == 3);
  CHECK(ids.size() == 3);
  CHECK(heights.size() == 3);
}

TEST_CASE("ViterbiPath dump_data_in_range subset") {
  ViterbiPath path(0);
  std::vector<int> hets = {};
  path.append(0, 1, 10.0, hets);
  path.append(3, 2, 20.0, hets);
  path.append(6, 3, 30.0, hets);

  std::vector<int> positions = {100, 200, 300, 400, 500, 600, 700, 800, 900};
  path.map_positions(positions);

  // Request range that covers only second segment
  auto [starts, ids, heights] = path.dump_data_in_range(400, 600);
  CHECK(starts.size() >= 1);
  // The first returned start should be 400
  CHECK(starts[0] == 400);
}

// === ViterbiState tests ===

TEST_CASE("ViterbiState construction") {
  std::vector<int> samples = {0, 1, 2};
  ViterbiState state(5, samples);
  CHECK(state.target_id == 5);
  CHECK(state.best_match == 0);
  CHECK(state.sites_processed == 0);
}

TEST_CASE("ViterbiState construction requires non-empty samples") {
  std::vector<int> empty_samples;
  CHECK_THROWS(ViterbiState(0, empty_samples));
}

TEST_CASE("ViterbiState process_site basic") {
  // 3 reference samples + target
  std::vector<int> samples = {0, 1, 2};
  ViterbiState state(3, samples);

  // Genotype vector: all samples + target
  // sample 0=0, sample 1=1, sample 2=0, target=1
  std::vector<int> geno = {0, 1, 0, 1};

  double rho = 5.0;    // recombination penalty
  double rho_c = 0.01; // non-recombination penalty
  double mu = 3.0;     // mutation penalty
  double mu_c = 0.01;  // non-mutation penalty

  state.process_site(geno, rho, rho_c, mu, mu_c);
  CHECK(state.sites_processed == 1);

  // Best match should be sample 1 (matches target allele)
  CHECK(state.best_match == 1);
}

TEST_CASE("ViterbiState process multiple sites and traceback") {
  std::vector<int> samples = {0, 1};
  ViterbiState state(2, samples);

  // Process several sites where sample 0 always matches target
  for (int i = 0; i < 5; i++) {
    std::vector<int> geno = {1, 0, 1}; // sample0=1, sample1=0, target=1
    state.process_site(geno, 5.0, 0.01, 3.0, 0.01);
  }

  CHECK(state.sites_processed == 5);

  auto path = state.traceback();
  CHECK(path.target_id == 2);
  CHECK(path.size() >= 1);
  // Best path should mostly copy sample 0
  CHECK(path.sample_ids[0] == 0);
}

TEST_CASE("ViterbiState prune reduces branch count") {
  std::vector<int> samples = {0, 1, 2, 3};
  ViterbiState state(4, samples);

  // Process enough sites to create branches
  for (int i = 0; i < 20; i++) {
    std::vector<int> geno;
    for (int s = 0; s < 5; s++) {
      geno.push_back(i % 2 == 0 ? s % 2 : (s + 1) % 2);
    }
    state.process_site(geno, 2.0, 0.5, 1.5, 0.1);
  }

  int branches_before = state.count_branches();
  state.prune();
  int branches_after = state.count_branches();

  // Prune should not increase branch count
  CHECK(branches_after <= branches_before);
  // Should still have at least as many branches as samples
  CHECK(branches_after >= static_cast<int>(samples.size()));
}

TEST_CASE("ViterbiState traceback produces valid path") {
  std::vector<int> samples = {0, 1};
  ViterbiState state(2, samples);

  // Alternating genotypes to force recombinations
  for (int i = 0; i < 10; i++) {
    std::vector<int> geno;
    if (i < 5) {
      geno = {1, 0, 1}; // target matches sample 0
    } else {
      geno = {0, 1, 1}; // target matches sample 1
    }
    state.process_site(geno, 1.0, 0.5, 2.0, 0.01);
  }

  auto path = state.traceback();
  CHECK(path.size() >= 1);
  // Segments should be ordered
  for (int i = 1; i < path.size(); i++) {
    CHECK(path.segment_starts[i] > path.segment_starts[i - 1]);
  }
}

TEST_CASE("ViterbiState set_samples updates candidate set") {
  std::vector<int> samples = {0, 1, 2, 3, 4};
  ViterbiState state(5, samples);

  // Process a few sites
  for (int i = 0; i < 3; i++) {
    std::vector<int> geno = {1, 0, 1, 0, 1, 1};
    state.process_site(geno, 2.0, 0.5, 1.5, 0.1);
  }

  // Reduce to subset
  std::unordered_set<int> new_samples = {0, 2};
  state.set_samples(new_samples);

  // Sample_ids should now contain new_samples + best_match
  CHECK(state.sample_ids.size() <= 4); // at most new_samples + best_match + margin
}
