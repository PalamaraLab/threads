This directory contains demographic models for the populations from the 1000 Genomes Project that may be used with Threads.

These demographies were inferred using [smc++](https://github.com/popgenmethods/smcpp) and shared with the [pyrho repository](https://github.com/popgenmethods/pyrho/blob/master/smcpp_popsizes_1kg.csv).

They have been altered such that time reflects generations, not years, and effective population sizes are measured in haploids instead of diploids.

This directory also contains a sample demographic model for a fixed Ne=10,000-sized population in `Ne10000.demo`.