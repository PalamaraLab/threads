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

#include "Matcher.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <unordered_set>
#include <vector>

namespace {

// Create evenly spaced genetic positions
std::vector<double> linear_positions(int n_sites, double start = 0.0, double step = 0.01) {
  std::vector<double> pos;
  pos.reserve(n_sites);
  for (int i = 0; i < n_sites; i++) {
    pos.push_back(start + i * step);
  }
  return pos;
}

} // namespace

TEST_CASE("Matcher construction basic") {
  int n_samples = 10;
  int n_sites = 100;
  auto positions = linear_positions(n_sites);

  Matcher m(n_samples, positions, 0.01, 0.5, 4, 2);
  CHECK(m.num_samples == n_samples);
  CHECK(m.num_sites == n_sites);
}

TEST_CASE("Matcher construction requires >= 3 sites") {
  auto pos2 = linear_positions(2);
  CHECK_THROWS_WITH(Matcher(5, pos2, 0.01, 0.5, 4, 2),
                     Catch::Matchers::ContainsSubstring("Need at least 3 sites"));
}

TEST_CASE("Matcher construction requires increasing positions") {
  std::vector<double> bad_pos = {0.0, 0.5, 0.3, 0.8};
  CHECK_THROWS_WITH(Matcher(5, bad_pos, 0.01, 0.5, 4, 2),
                     Catch::Matchers::ContainsSubstring("strictly increasing"));
}

TEST_CASE("Matcher process_site with binary genotypes") {
  int n_samples = 20;
  int n_sites = 50;
  auto positions = linear_positions(n_sites, 0.0, 0.02);

  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  // Process all sites with alternating genotypes
  for (int site = 0; site < n_sites; site++) {
    std::vector<int> geno(n_samples);
    for (int s = 0; s < n_samples; s++) {
      geno[s] = (s + site) % 2;
    }
    m.process_site(geno);
  }

  auto matches = m.get_matches();
  CHECK(matches.size() > 0);
}

TEST_CASE("Matcher rejects wrong genotype size") {
  int n_samples = 10;
  auto positions = linear_positions(20, 0.0, 0.02);
  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  // Wrong size genotype
  std::vector<int> bad_geno(5, 0);
  CHECK_THROWS_WITH(m.process_site(bad_geno),
                     Catch::Matchers::ContainsSubstring("invalid genotype vector size"));
}

TEST_CASE("Matcher rejects invalid alleles") {
  int n_samples = 5;
  auto positions = linear_positions(10, 0.0, 0.02);
  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  std::vector<int> bad_geno = {0, 1, 0, 2, 0}; // 2 is invalid
  CHECK_THROWS_WITH(m.process_site(bad_geno),
                     Catch::Matchers::ContainsSubstring("invalid genotype"));
}

TEST_CASE("Matcher process_site rejects extra sites") {
  int n_samples = 5;
  auto positions = linear_positions(4, 0.0, 0.02);
  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  std::vector<int> geno(5, 0);
  for (int i = 0; i < 4; i++) {
    m.process_site(geno);
  }
  CHECK_THROWS_WITH(m.process_site(geno),
                     Catch::Matchers::ContainsSubstring("all sites have already been processed"));
}

TEST_CASE("Matcher sorting is a valid permutation") {
  int n_samples = 10;
  int n_sites = 30;
  auto positions = linear_positions(n_sites, 0.0, 0.02);

  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  for (int site = 0; site < n_sites; site++) {
    std::vector<int> geno(n_samples);
    for (int s = 0; s < n_samples; s++) {
      geno[s] = (s * 3 + site) % 2;
    }
    m.process_site(geno);
  }

  auto sorting = m.get_sorting();
  CHECK(static_cast<int>(sorting.size()) == n_samples);

  // Check it's a valid permutation
  std::unordered_set<int> seen;
  for (int v : sorting) {
    CHECK(v >= 0);
    CHECK(v < n_samples);
    seen.insert(v);
  }
  CHECK(static_cast<int>(seen.size()) == n_samples);
}

TEST_CASE("Matcher cm_positions") {
  int n_samples = 10;
  int n_sites = 100;
  auto positions = linear_positions(n_sites, 0.0, 0.01);

  Matcher m(n_samples, positions, 0.01, 0.5, 4, 1);

  // Process all sites
  for (int site = 0; site < n_sites; site++) {
    std::vector<int> geno(n_samples);
    for (int s = 0; s < n_samples; s++) {
      geno[s] = s % 2;
    }
    m.process_site(geno);
  }

  auto cms = m.cm_positions();
  CHECK(cms.size() > 0);
  // Should be non-decreasing
  for (std::size_t i = 1; i < cms.size(); i++) {
    CHECK(cms[i] >= cms[i - 1]);
  }
}

TEST_CASE("MatchGroup construction") {
  MatchGroup mg(10, 0.5);
  CHECK(mg.num_samples == 10);
  CHECK(mg.cm_position == 0.5);
  CHECK(mg.match_candidates_counts.size() == 10);
}

TEST_CASE("MatchGroup from targets and matches") {
  std::vector<int> targets = {0, 1, 2};
  std::vector<std::unordered_set<int>> matches = {{}, {0}, {0, 1}};
  MatchGroup mg(targets, matches, 1.0);

  CHECK(mg.match_candidates.size() == 3);
  CHECK(mg.match_candidates.at(0).size() == 0);
  CHECK(mg.match_candidates.at(1).size() == 1);
  CHECK(mg.match_candidates.at(2).size() == 2);
}

TEST_CASE("MatchGroup clear") {
  MatchGroup mg(5, 0.0);
  mg.clear();
  CHECK(mg.match_candidates.empty());
  CHECK(mg.match_candidates_counts.empty());
  CHECK(mg.top_four_maps.empty());
}
