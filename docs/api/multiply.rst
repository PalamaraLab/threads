Genotype-Matrix Multiplication
==============================

Multiply the implicit genotype matrix by vectors without materializing the full
matrix. Three strategies are available, from simplest to fastest.

Dense Multiply (baseline)
-------------------------

Materializes the full genotype matrix from threading instructions, then computes
standard matrix-vector products. O(n * m) per call and O(n * m) memory.
Useful as a correctness reference and for small problems.

.. method:: ThreadingInstructions.left_multiply(x, diploid=False, normalize=False)

   Compute G.T @ x (dense). If ``diploid``, uses diploid dosage matrix
   (pairs haploid samples). If ``normalize``, uses mean-centered
   variance-normalized matrix.

.. method:: ThreadingInstructions.right_multiply(x, diploid=False, normalize=False)

   Compute G @ x (dense). Same options as ``left_multiply``.

RLE Multiply (sparse baseline)
------------------------------

Stores each sample's genotype as runs of consecutive 1s — equivalent to a sparse
matrix representation. Simpler than DAG but linear in m.

**Complexity:** O(m + total_1_runs) per multiply.

.. method:: ThreadingInstructions.prepare_rle_multiply()

   Precompute 1-run (start, end) pairs per sample. Must be called before RLE
   multiply methods.

.. method:: ThreadingInstructions.right_multiply_rle(x, diploid=False, normalize=False)

   Compute G @ x via RLE. Accepts 1D ``(num_sites,)`` or 2D ``(num_sites, k)``
   input. Returns numpy array of matching shape.

.. method:: ThreadingInstructions.left_multiply_rle(x, diploid=False, normalize=False)

   Compute G.T @ x via RLE. Accepts 1D ``(num_samples,)`` or 2D
   ``(num_samples, k)`` input. Returns numpy array of matching shape.

DAG Multiply (fast, sublinear)
------------------------------

Preprocesses threading instructions into a compact directed acyclic graph of
mutation sets, following the ARG multiplication approach used in
arg-needle-lib (ARGMatMult). Segments without mismatches
are collapsed, yielding a DAG whose size grows sublinearly with n (empirically
~n^0.64). The genomic region is split into independent chunks for OpenMP
parallelism.

**Complexity:** O(|DAG edges| + n) per multiply — sublinear when the DAG is
small relative to n * m.

.. method:: ThreadingInstructions.prepare_dag_multiply(num_chunks=0)

   Build the chunked mutation-set DAG. ``num_chunks`` defaults to
   ``OMP_NUM_THREADS``. Must be called before DAG multiply methods.

.. method:: ThreadingInstructions.right_multiply_dag(x, diploid=False, normalize=False)

   Compute G @ x via DAG. Accepts 1D ``(num_sites,)`` input.
   Returns numpy array ``(num_samples,)``.

.. method:: ThreadingInstructions.left_multiply_dag(x, diploid=False, normalize=False)

   Compute G.T @ x via DAG. Accepts 1D ``(num_samples,)`` input.
   Returns numpy array ``(num_sites,)``.

.. method:: ThreadingInstructions.right_multiply_dag_batch(X, k)

   Batch right multiply: X is ``(num_sites, k)`` flat array.
   Returns ``(num_samples, k)``. Columns distributed across OMP threads.

.. method:: ThreadingInstructions.left_multiply_dag_batch(X, k)

   Batch left multiply: X is ``(num_samples, k)`` flat array.
   Returns ``(num_sites, k)``. Columns distributed across OMP threads.

Diploid and Normalize Flags
----------------------------

All multiply methods accept ``diploid`` and ``normalize`` keyword arguments:

- **diploid**: Operates on diploid dosage (sum of haploid pairs) instead of
  haploid genotypes. Right multiply reduces input from n to n/2 diploid
  individuals; left multiply expands output from n/2 back to n sites.

- **normalize**: Mean-centers and variance-normalizes the genotype matrix.
  Uses empirical dosage variance for diploid mode. Equivalent to multiplying
  by ``(G - mean) / std``.

These flags are applied as pre/post-processing wrappers around the raw multiply,
matching the ARGMatMult convention.
