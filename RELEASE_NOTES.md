# threads Release Notes

## [Unreleased]

### Added

- Add left_multiplication and right_multiplication to ThreadingInstructions (#99)
- OpenMP parallelism for Li-Stephens Viterbi and hets phases across targets (#127)
- NumPy batch API for zero-copy genotype transfer from Python to C++ (#127)

### Changed

- Build wheels on macOS 14 for arm64 and macOS 15 for x86_64 (#108)
- Replace hash map traceback storage with deque arena in Viterbi (#127)
- Replace hash map state lookup with flat vectors in ThreadsLowMem (#127)
- Boolean array neighbor search in Matcher replaces red-black tree (#127)

### Fixed

- Memory leak: raw HMM pointer replaced with unique_ptr (#127)
- Sign-extension in pair_key/coord_id_key hash functions (#127)
- Bounds check after array access in ThreadsFastLS (#127)
- exit(1) calls replaced with exceptions for proper Python error propagation (#127)

## [0.2.1] - 2025-06-03

### Added

- VCF conversion backwards compatibility (#76)
- Add Linux ARM platform to Python wheels (#82)

## [0.2.0] - 2025-04-29

### Added

- Optimised imputation (#28)
- Variant mapping (#38)
- Upgrade to numpy 2.0 (#39)
- Consistent SHAPEIT format for regions (#42)
- Ability to write impute output directly to stdout (#42)
- Allele age estimation and data consistency (#51)

### Fixed

- Fix os.sched_getaffinity macOS error (#49)

## [0.1.0] - 2024-07-04

### Added

- Initial version of threads
