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
//
// Benchmarks for core algorithms. Run with:
//   ./unit_tests "[benchmark]" --benchmark-samples 5
// or without Catch2 benchmark (we do manual timing):
//   ./unit_tests "[benchmark]"

#include "Demography.hpp"
#include "HMM.hpp"
#include "Matcher.hpp"
#include "ViterbiLowMem.hpp"

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#ifdef __APPLE__
#include <mach/mach.h>
static size_t get_resident_memory_bytes() {
  struct mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS) {
    return info.resident_size;
  }
  return 0;
}
#elif defined(__linux__)
#include <fstream>
#include <string>
#include <unistd.h>
static size_t get_resident_memory_bytes() {
  std::ifstream f("/proc/self/statm");
  size_t pages;
  f >> pages; // total
  f >> pages; // resident
  return pages * sysconf(_SC_PAGESIZE);
}
#else
static size_t get_resident_memory_bytes() { return 0; }
#endif

namespace {

struct BenchResult {
  double elapsed_ms;
  size_t mem_before;
  size_t mem_after;

  double mem_delta_mb() const {
    return static_cast<double>(mem_after - mem_before) / (1024.0 * 1024.0);
  }

  void print(const char* label) const {
    std::cout << "  [BENCH] " << label << ": " << std::fixed << std::setprecision(2) << elapsed_ms
              << " ms, mem delta: " << std::setprecision(2) << mem_delta_mb() << " MB" << std::endl;
  }
};

// Generate deterministic genotypes
std::vector<std::vector<int>> generate_genotypes(int n_samples, int n_sites, unsigned seed = 42) {
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(0, 1);
  std::vector<std::vector<int>> genos(n_sites, std::vector<int>(n_samples));
  for (int s = 0; s < n_sites; s++) {
    for (int i = 0; i < n_samples; i++) {
      genos[s][i] = dist(rng);
    }
  }
  return genos;
}

std::vector<double> linear_positions(int n_sites, double start, double step) {
  std::vector<double> pos;
  pos.reserve(n_sites);
  for (int i = 0; i < n_sites; i++) {
    pos.push_back(start + i * step);
  }
  return pos;
}

} // namespace

TEST_CASE("Benchmark: HMM construction", "[benchmark]") {
  const int n_sites = 10000;
  const int K = 64;
  Demography demo({10000.0}, {0.0});
  auto bp = std::vector<double>(n_sites, 100.0);
  auto cm = std::vector<double>(n_sites, 0.001);

  size_t mem_before = get_resident_memory_bytes();
  auto t0 = std::chrono::high_resolution_clock::now();

  HMM hmm(demo, bp, cm, 1.4e-8, K);

  auto t1 = std::chrono::high_resolution_clock::now();
  size_t mem_after = get_resident_memory_bytes();

  BenchResult r{
      std::chrono::duration<double, std::milli>(t1 - t0).count(), mem_before, mem_after};
  r.print("HMM construction (10K sites, K=64)");

  // Sanity
  CHECK(hmm.trellis.size() == n_sites);
}

TEST_CASE("Benchmark: HMM breakpoints", "[benchmark]") {
  const int n_sites = 5000;
  const int K = 64;
  Demography demo({10000.0}, {0.0});
  auto bp = std::vector<double>(n_sites, 100.0);
  auto cm = std::vector<double>(n_sites, 0.001);
  HMM hmm(demo, bp, cm, 1.4e-8, K);

  // Mixed observation pattern
  std::mt19937 rng(123);
  std::uniform_int_distribution<int> dist(0, 4);
  std::vector<bool> obs(n_sites);
  for (int i = 0; i < n_sites; i++) {
    obs[i] = dist(rng) == 0; // ~20% het rate
  }

  auto t0 = std::chrono::high_resolution_clock::now();

  auto bps = hmm.breakpoints(obs, 0);

  auto t1 = std::chrono::high_resolution_clock::now();
  BenchResult r{std::chrono::duration<double, std::milli>(t1 - t0).count(), 0, 0};
  r.print("HMM breakpoints (5K sites, K=64)");

  CHECK(bps.size() >= 1);
}

TEST_CASE("Benchmark: Matcher process_site", "[benchmark]") {
  const int n_samples = 1000;
  const int n_sites = 500;
  auto positions = linear_positions(n_sites, 0.0, 0.02);
  auto genos = generate_genotypes(n_samples, n_sites);

  size_t mem_before = get_resident_memory_bytes();
  auto t0 = std::chrono::high_resolution_clock::now();

  Matcher m(n_samples, positions, 0.01, 0.5, 4, 2);
  for (int s = 0; s < n_sites; s++) {
    m.process_site(genos[s]);
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  size_t mem_after = get_resident_memory_bytes();

  BenchResult r{
      std::chrono::duration<double, std::milli>(t1 - t0).count(), mem_before, mem_after};
  r.print("Matcher process_site (1K samples, 500 sites)");

  CHECK(m.get_sorting().size() == n_samples);
}

TEST_CASE("Benchmark: ViterbiState process_site", "[benchmark]") {
  const int n_ref = 100;
  const int target_id = n_ref;
  const int n_sites = 2000;
  const int n_samples_total = n_ref + 1;

  std::vector<int> ref_samples;
  for (int i = 0; i < n_ref; i++) {
    ref_samples.push_back(i);
  }

  auto genos = generate_genotypes(n_samples_total, n_sites);

  size_t mem_before = get_resident_memory_bytes();
  auto t0 = std::chrono::high_resolution_clock::now();

  ViterbiState state(target_id, ref_samples);
  for (int s = 0; s < n_sites; s++) {
    state.process_site(genos[s], 3.0, 0.01, 2.0, 0.01);
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  size_t mem_after = get_resident_memory_bytes();

  BenchResult r{
      std::chrono::duration<double, std::milli>(t1 - t0).count(), mem_before, mem_after};
  r.print("ViterbiState process_site (100 refs, 2K sites)");

  CHECK(state.sites_processed == n_sites);
}

TEST_CASE("Benchmark: ViterbiState prune", "[benchmark]") {
  const int n_ref = 50;
  const int target_id = n_ref;
  const int n_sites = 1000;
  const int n_samples_total = n_ref + 1;

  std::vector<int> ref_samples;
  for (int i = 0; i < n_ref; i++) {
    ref_samples.push_back(i);
  }

  auto genos = generate_genotypes(n_samples_total, n_sites);

  ViterbiState state(target_id, ref_samples);
  for (int s = 0; s < n_sites; s++) {
    state.process_site(genos[s], 3.0, 0.01, 2.0, 0.01);
  }

  int branches_before = state.count_branches();

  auto t0 = std::chrono::high_resolution_clock::now();
  state.prune();
  auto t1 = std::chrono::high_resolution_clock::now();

  int branches_after = state.count_branches();

  BenchResult r{std::chrono::duration<double, std::milli>(t1 - t0).count(), 0, 0};
  r.print("ViterbiState prune (50 refs, 1K sites)");

  std::cout << "    Branches: " << branches_before << " -> " << branches_after << std::endl;
  CHECK(branches_after <= branches_before);
}

TEST_CASE("Benchmark: ViterbiState traceback", "[benchmark]") {
  const int n_ref = 50;
  const int target_id = n_ref;
  const int n_sites = 1000;
  const int n_samples_total = n_ref + 1;

  std::vector<int> ref_samples;
  for (int i = 0; i < n_ref; i++) {
    ref_samples.push_back(i);
  }

  auto genos = generate_genotypes(n_samples_total, n_sites);

  ViterbiState state(target_id, ref_samples);
  for (int s = 0; s < n_sites; s++) {
    state.process_site(genos[s], 3.0, 0.01, 2.0, 0.01);
  }

  auto t0 = std::chrono::high_resolution_clock::now();
  auto path = state.traceback();
  auto t1 = std::chrono::high_resolution_clock::now();

  BenchResult r{std::chrono::duration<double, std::milli>(t1 - t0).count(), 0, 0};
  r.print("ViterbiState traceback (50 refs, 1K sites)");

  CHECK(path.size() >= 1);
  std::cout << "    Path segments: " << path.size() << std::endl;
}

TEST_CASE("Benchmark: Demography std_to_gen (many calls)", "[benchmark]") {
  Demography d({5000.0, 10000.0, 20000.0}, {0.0, 100.0, 500.0});

  const int n_calls = 1000000;
  double sum = 0.0;

  auto t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < n_calls; i++) {
    double t = static_cast<double>(i) / n_calls * 5.0;
    sum += d.std_to_gen(t);
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  BenchResult r{std::chrono::duration<double, std::milli>(t1 - t0).count(), 0, 0};
  r.print("Demography std_to_gen (1M calls)");

  CHECK(sum > 0.0); // prevent optimization
}
