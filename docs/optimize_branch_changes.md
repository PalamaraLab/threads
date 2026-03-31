# Optimize Branch: Detailed Change Description

**PR #5**: `optimize` → `main` on `threads-dev`
**Scope**: 43 files changed, +6,465 / -899 lines across 10 commits

---

## 1. Genotype-Matrix Multiply (`ThreadingInstructions`)

Three multiply strategies operate on the implicit genotype matrix encoded by threading instructions.

### 1.1 Dense multiply (baseline)

Materializes the full n×m genotype matrix, then computes standard matrix-vector products. O(n*m) time and memory. Useful as a correctness reference.

Caching methods avoid repeated `GenotypeIterator` traversals:
- `materialize_genotypes()` — caches the full integer genotype matrix.
- `materialize_normalized_haploid()` — caches `(g - mu) / std` as doubles.
- `materialize_normalized_diploid()` — caches diploid (paired-sum) standardized matrix.

### 1.2 RLE multiply (sparse baseline)

Stores each sample's genotype as runs of consecutive 1s. O(m + total_1_runs) per multiply. Simpler than DAG but scales linearly in m.

### 1.3 DAG multiply (fast, sublinear)

Preprocesses threading instructions into a compact directed acyclic graph of mutation sets, following the ARG multiplication approach used in arg-needle-lib (ARGMatMult). Segments without mismatches are collapsed (pass-through), yielding a DAG whose size grows sublinearly with n.

**`prepare_dag_multiply(num_chunks)`** (one-time precomputation):
- Builds a reference genome by tracing sample 0's threading chain.
- Splits the genomic region into `num_chunks` independent sub-DAGs (default: OMP_NUM_THREADS).
- For each chunk: identifies unique mutation sets via a reverse pass through threading chains, then builds parent→child and child→parent edge lists, mutation attachments, and individual→set mappings.

**`right_multiply_dag(x)`** — computes Gx:
- Bottom-up propagation through DAG sets: each set accumulates from its descendants, applies mutation corrections, then scatters to connected individuals.
- Per-chunk processing is embarrassingly parallel across OMP threads.

**`left_multiply_dag(x)`** — computes G^T x:
- Gathers individual contributions per set, propagates through the DAG, then scatters mutation-weighted sums to output sites.

**Batch variants** (`right_multiply_dag_batch`, `left_multiply_dag_batch`):
- Process k vectors via column-parallel OMP. DAG structure is shared read-only; per-thread data is private.

### 1.4 Diploid and normalize flags

All multiply methods accept `diploid` and `normalize` keyword arguments, applied as pre/post-processing wrappers (matching the ARGMatMult convention):
- **diploid**: right multiply reduces haploid pairs to dosage; left multiply expands back.
- **normalize**: mean-centers and variance-normalizes using empirical dosage variance.

### 1.5 Complexity

| Operation | Time | Memory |
|-----------|------|--------|
| Dense `left/right_multiply` | O(n * m) | O(n * m) cached matrix |
| RLE `left/right_multiply_rle` | O(m + total_1_runs) | O(total_runs) |
| DAG `prepare_dag_multiply` | O(n * m) one-time | O(\|DAG\|) |
| DAG `left/right_multiply_dag` | O(\|DAG edges\| + n) | O(n) per call |

---

## 2. OpenMP Parallelism (`ThreadsLowMem`)

The Li-Stephens Viterbi decoding and hets-counting phases are parallelized across target samples using OpenMP. This is the main inference bottleneck at scale.

### 2.1 What changed

Previously, `ThreadsLowMem` processed targets sequentially:
```cpp
// Old: sequential
for (const auto& genotype : genotypes) {
    process_site_viterbi_raw(genotype.data());
}
```

Now, batch methods restructure the loop to parallelize across targets:
```cpp
// New: parallel across targets
#pragma omp parallel for schedule(dynamic, 1)
for (int idx = 0; idx < n_active; ++idx) {
    // Each thread processes ALL sites for one target
    for (int s = 0; s < n_sites_batch; s++) {
        // ... Viterbi update for this target at this site
    }
}
```

Key design choice: the outer loop is over targets (not sites), because each target's Viterbi state is independent. This avoids synchronization between sites and gives each thread a contiguous block of work.

### 2.2 Affected methods

| Method | Parallelism |
|--------|-------------|
| `process_all_sites_viterbi()` | `#pragma omp parallel for` over targets |
| `process_all_sites_viterbi_flat()` | Same, with flat `int32_t*` input |
| `process_all_sites_hets()` | `#pragma omp parallel for` over targets |
| `process_all_sites_hets_flat()` | Same, with flat input |
| `traceback()` | `#pragma omp parallel for` over targets |
| `date_segments()` | `#pragma omp parallel for` over targets |
| `prune()` | `#pragma omp parallel for` over HMM instances |

### 2.3 Precomputed per-site constants

To avoid repeated floating-point computation inside the parallel loop, `k_per_site` and `l_per_site` vectors are precomputed once during initialization:
```cpp
k_per_site[i] = 2.0 * 0.01 * cm_sizes[i];   // recombination
l_per_site[i] = 2.0 * mutation_rate * bp_sizes[i];  // mutation
```

### 2.4 Build system

`CMakeLists.txt` adds optional OpenMP detection:
```cmake
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(threads_arg_python_bindings PRIVATE OpenMP::OpenMP_CXX)
endif()
```
The code compiles and runs correctly without OpenMP (the `#pragma` is ignored).

---

## 3. Matcher Optimizations

### 3.1 Boolean array neighbor search

The PBWT matching step previously used `std::set` (red-black tree) to find nearest neighbors in the permutation order. This was O(log n) per insertion and query.

Replaced with a `std::vector<char>` boolean array + sequential scan:
```cpp
std::vector<char> inserted(num_samples, 0);
// ...
while (left >= 0 && !inserted[left]) left--;  // scan left
while (right < num_samples && !inserted[right]) right++;  // scan right
```

This is O(gap_size) per query but much faster in practice due to cache locality and branch prediction on the boolean array. The neighborhood_size is small (typically 4), so only a few neighbors are needed.

### 3.2 Raw pointer API

Added `process_site_raw(const int* genotype)` which takes a raw pointer instead of `std::vector<int>`. This eliminates a vector copy on every site when called from the NumPy batch API. The existing `process_site(const std::vector<int>&)` now delegates to `process_site_raw`.

### 3.3 NumPy batch API

`process_all_sites_flat(const int32_t* data, int n_sites, int n_haps)` accepts a flat NumPy array from Python (zero-copy via pybind11) and processes all sites without per-site Python→C++ round trips. On the Python side:
```python
if HAS_NUMPY_API:
    matcher.process_all_sites_numpy(all_genotypes[ac_mask])
else:
    for g in all_genotypes[ac_mask]:
        matcher.process_site(g)
```

---

## 4. Bug Fixes

### 4.1 `exit(1)` → exceptions

20 instances of `exit(1)` in the C++ codebase were replaced with `throw std::runtime_error(...)`. The old `exit(1)` calls would terminate the Python process without any opportunity for cleanup or error handling. Now errors propagate as Python exceptions via pybind11.

### 4.2 Memory leak in HMM allocation

`ThreadsLowMem` allocated `ViterbiState` objects with `new` and stored raw pointers:
```cpp
// Old
HMM* hmm = new HMM(...);
hmm_ptrs.push_back(hmm);  // never deleted
```
Replaced with `std::unique_ptr`:
```cpp
// New
hmm_vec.push_back(std::make_unique<ViterbiState>(...));
hmm_ptrs.push_back(hmm_vec.back().get());
```

### 4.3 Sign-extension in hash functions

`pair_key` combined two `int` values into a `long long` using bit shifts, but the upper int could sign-extend into the high bits, causing hash collisions:
```cpp
// Old (buggy)
long long pair_key(int a, int b) { return ((long long)a << 32) | b; }
// New
long long pair_key(int a, int b) { return ((long long)a << 32) | (uint32_t)b; }
```

### 4.4 Self-referencing target bugs

In `AlleleAges.cpp`: carrier chain traversal and trace path could enter infinite loops when a sample's target was itself (self-referencing segments). Added guards:
```cpp
if (target == sample) break;  // self-ref: stop traversal
```

In `DataConsistency.cpp`: self-referencing targets caused silent mismatch skipping. Fixed by substituting sample 0 when a self-ref is encountered.

### 4.5 Sign bug in rho computation

Visible at scale — the rho (recombination) parameter had a sign error that produced incorrect Viterbi paths for large sample sizes. Fixed in the precomputed `k_per_site` initialization.

---

## 5. Dependency Cleanup

### 5.1 Removed from core dependencies

| Package | Reason removed | Replacement |
|---------|---------------|-------------|
| `ray` | Heavy (~1s import), required distributed scheduler | `multiprocessing` (stdlib) |
| `click` | ~224ms import time for CLI | `argparse` (stdlib) |
| `pandas` | Used only for reading tab-delimited files | Direct file parsing with stdlib |
| `numba` | JIT compilation of forward/backward algorithm | Ported to C++ (`ForwardBackward.cpp`) |
| `xarray` | Unused in current code | Removed |
| `h5py` | Used for serialization of ThreadingInstructions | Ported to C++ (`ThreadsIO.cpp`) |
| `tszip` | Only needed for `.tsz` conversion output | Moved to optional `[convert]` group |

### 5.2 Optional dependency groups

```toml
[project.optional-dependencies]
dev = ["pytest", "h5py"]
convert = ["tszip"]
impute = ["pandas", "scipy"]
normalize = ["msprime"]
all = ["tszip", "msprime", "pandas", "scipy"]
```

### 5.3 Lazy imports with error messages

All optional imports are now lazy with descriptive error messages:
```python
try:
    from threads_arg.serialization import serialize_instructions
except ImportError:
    # Only imported when --out is specified
    pass
```
This means `threads infer` can run without `h5py`, `tszip`, or `msprime` installed unless those specific features are requested.

---

## 6. Inference Pipeline Changes (`infer.py`)

### 6.1 Read-once genotype loading

Previously, genotypes were read from the pgen file multiple times (once for singleton filtering, once for matching, once for Viterbi, once for hets). Now:
```python
all_genotypes = read_all_genotypes(pgen)  # single read into numpy array
```
All subsequent stages index into this array.

### 6.2 Parallelism model change

**Old**: `ray` distributed the Li-Stephens Viterbi across processes, with each process handling a batch of targets across all sites. Required serializing match data to each worker.

**New**: Two paths depending on build:
1. **OpenMP path** (preferred): Single process, all genotypes in memory. C++ `process_all_sites_viterbi()` handles parallelism across targets via OpenMP. No serialization overhead.
2. **Multiprocessing fallback**: When OpenMP is unavailable, uses `multiprocessing.Pool` with the same batch-per-target strategy as the old ray implementation, but without the ray dependency.

### 6.3 Consistency check moved to C++

The data consistency check (verifying inferred threading reproduces input genotypes) was previously a Python loop calling `ConsistencyWrapper` per-site. Now delegated to `run_consistency()` which is a single C++ call processing all sites internally.

---

## 7. Python-Side Optimizations

### 7.1 CLI: `click` → `argparse`

The entire CLI (`__main__.py`) was rewritten from `click` decorators to `argparse` subparsers. Saves ~224ms of import time on every invocation.

### 7.2 File I/O: `pandas` → stdlib

`read_map_file()`, `_read_pgen_physical_positions()`, and `read_variant_metadata()` in `utils.py` were rewritten to parse files directly with `open()` + `split()` instead of `pd.read_table()`. A lightweight `VariantMetadata` class replaces `pandas.DataFrame` for variant metadata.

### 7.3 Forward/backward algorithm: `numba` → C++

The forward/backward Li-Stephens algorithm in `fwbw.py` was JIT-compiled with numba. Now ported to C++ in `ForwardBackward.cpp` (95 lines) and exposed via pybind11. This eliminates the numba dependency and its ~2s JIT compilation on first call.

### 7.4 Serialization: `h5py` → C++

HDF5 reading/writing of ThreadingInstructions moved to `ThreadsIO.cpp` (338 lines). The Python `serialization.py` now calls C++ `load_instructions()` / `save_instructions()` when available, falling back to the Python h5py implementation if needed.

### 7.5 Allele ages

- Removed unnecessary `np.array()` conversion in per-site loop (was creating a new array from an existing array on every iteration).
- Replaced `multiprocessing.Manager().dict()` (which serializes every dict update through a manager process) with `pool.imap()` return values.
- Added `target < 0` guard for out-of-bounds crash when carrier chain encounters invalid targets.

### 7.6 Impute

- Vectorized emission probability computation in `fwbw.py`.
- Optimized sparse posterior construction in `impute.py` using `scipy.sparse.csr_array`.
- Improved VCF writing buffering.

---

## 8. New C++ Modules

### 8.1 `ForwardBackward.cpp` / `.hpp` (95 + 42 lines)

C++ port of the forward/backward Li-Stephens HMM, replacing the numba-JIT Python implementation. Exposed via pybind11 as `forwards_ls_hap_cpp` and `backwards_ls_hap_cpp`.

### 8.2 `ThreadsIO.cpp` / `.hpp` (338 + 21 lines)

C++ HDF5 reader/writer for `ThreadingInstructions` serialization. Uses the HDF5 C API directly. Replaces the Python `h5py`-based serialization. Exposed as `load_instructions()` and `save_instructions()`.

---

## 9. Python Bindings (`threads_arg_pybind.cpp`)

New bindings added:

| Binding | Description |
|---------|-------------|
| `prepare_dag_multiply(num_chunks)` | Precompute chunked DAG for sublinear multiply |
| `right_multiply_dag(x, diploid, normalize)` | DAG Gx, with GIL release |
| `left_multiply_dag(x, diploid, normalize)` | DAG G^T x, with GIL release |
| `prepare_rle_multiply()` | Precompute RLE sparse baseline |
| `right_multiply_rle(x, diploid, normalize)` | RLE Gx, with GIL release |
| `left_multiply_rle(x, diploid, normalize)` | RLE G^T x, with GIL release |
| `left_multiply(x, diploid, normalize)` | Dense baseline G^T x |
| `right_multiply(x, diploid, normalize)` | Dense baseline Gx |
| `materialize_genotypes()` | Cache dense genotype matrix |
| `materialize_normalized_haploid()` | Cache standardized haploid matrix |
| `materialize_normalized_diploid()` | Cache standardized diploid matrix |
| `process_all_sites_numpy(arr)` | Matcher: batch site processing from NumPy |
| `run_consistency(...)` | Bulk data consistency check |
| `load_instructions(path)` | C++ HDF5 deserialization |
| `save_instructions(ti, path)` | C++ HDF5 serialization |

All multiply and Viterbi bindings use `py::call_guard<py::gil_scoped_release>()`, enabling Python-level thread parallelism (e.g., `ThreadPoolExecutor`) to dispatch multiple C++ calls concurrently.

---

## 10. Testing

### 10.1 C++ tests

69 new Catch2 unit tests (444 assertions) in `test_benchmark.cpp` and `test_allele_ages.cpp`, covering:
- Tree multiply correctness (vs dense multiply baseline)
- Batch multiply consistency
- Range-restricted multiply
- Memory benchmarks (resident memory tracking)
- Allele age estimation edge cases

### 10.2 Python tests (~2,500 lines)

| Test file | Coverage |
|-----------|----------|
| `test_infer.py` | Utility functions, Matcher, ThreadsLowMem, TI construction, serialization round-trips |
| `test_allele_ages.py` | AgeEstimator, ConsistencyWrapper, synthetic TI properties |
| `test_impute_correctness.py` | Forward/backward algorithm, emission probabilities, sparse posteriors, VCF output |
| `test_convert.py` | .threads → .argn / .tsz conversion |
| `test_map.py` | Mutation mapping to ARG |
| `test_normalization.py` | Genotype normalization pipeline |
| `test_vcf.py` | VCF writing, genotype formatting, metadata handling |

### 10.3 Benchmarks

| Benchmark | Purpose |
|-----------|---------|
| `bench_compression_multiply.py` | Compare threads vs GRG vs RePair multiply speed and compression |
| `bench_impute.py` | End-to-end imputation timing |
| `bench_impute_scaling.py` | Imputation scaling with panel size |
| `build_cache.py` | Build cached simulation data for benchmarks |

---

## 11. Removed Code

- **`phase.py`**: Dead module with broken imports and no CLI registration. Removed entirely.
- **`compare_vs_genrepair.py`**: Superseded by `bench_paper.py` in threads-tools.
- **Ray dependency and distributed worker code** in `infer.py`.
