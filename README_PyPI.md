# Threads

This software implements the Threads algorithm, described in

`Á. F. Gunnarsson, J. Zhu, B. C. Zhang, Z. Tsangalidou, A. Allmont, P. Palamara. A scalable approach for genome-wide inference of ancestral recombination graphs. bioRxiv, 2024.`

The user manual for threads can be found [here](https://palamaralab.github.io/software/threads/).

## Installation

Prebuilt CPython wheels are available for Linux (compatible with glibc ≥ 2.28) and macOS (built on macOS 15 for x86_64 and macOS 14 for arm64).

| Platform \ CPython          | ≤3.8 | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|-----------------------------| ---- | --- | ---- | ---- | ---- | ---- | ---- |
| Linux x86_64                | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| Linux aarch64               | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Intel (x86_64)        | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Apple Silicon (arm64) | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |

```sh
pip install threads_arg
```
