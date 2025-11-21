# threads Release Notes

## [Unreleased]

### Changed

- Build wheels on macOS 14 for arm64 and macOS 15 for x86_64 (#108)

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
