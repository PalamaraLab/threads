Genotype-Matrix Multiplication
==============================

Multiply the implicit genotype matrix by vectors without materializing the full
matrix. Two strategies are available: tree shuttle (preferred) and RLE.

Tree Shuttle
------------

Propagates values through the interval-tree structure of threading instructions.
Both directions use O(n) per-call memory and OpenMP multithreading.

**Right multiply:** O((n * n_intervals + mismatches) / T) per call, where T is
the number of threads. Region-chunked parallelism with O(n*T) memory.

**Left multiply:** O((m + mismatches + segments * depth) / T) per call.
Uses an incremental O(n) subtree weight vector (not a weight matrix), achieving
sublinear scaling in n. Empirically ~n^0.53.

.. method:: ThreadingInstructions.prepare_tree_multiply()

   Precompute tree-shuttle structures: interval decomposition, mismatch signs
   via lazy chain tracing, inverted mismatch index, interval change list, and
   precomputed ancestor chains. O(n * S * d) time where S = segments per sample
   and d = tree depth. No genotype matrix is materialized.

.. method:: ThreadingInstructions.right_multiply_tree(X)

   Compute G @ X. Accepts 1D array/list ``(num_sites,)`` or 2D array
   ``(num_sites, k)``. Returns numpy array of matching shape:
   ``(num_samples,)`` or ``(num_samples, k)``.

.. method:: ThreadingInstructions.left_multiply_tree(X)

   Compute G.T @ X. Accepts 1D array/list ``(num_samples,)`` or 2D array
   ``(num_samples, k)``. Returns numpy array of matching shape:
   ``(num_sites,)`` or ``(num_sites, k)``.

.. method:: ThreadingInstructions.right_multiply_tree_range(x, site_start, site_end)

   Compute G[site_start:site_end, :] @ x. Returns list of length ``num_samples``.

.. method:: ThreadingInstructions.left_multiply_tree_range(x, site_start, site_end)

   Compute G[site_start:site_end, :].T @ x. Returns list of length
   ``site_end - site_start``.

RLE (Run-Length Encoding)
-------------------------

Stores each sample's genotype as runs of 1s---equivalent to a sparse matrix
representation. Simpler than tree shuttle but dominated by the O(m) prefix-sum
term.

**Complexity:** O(m + total_1_runs) per multiply.

.. method:: ThreadingInstructions.prepare_rle_multiply()

   Precompute 1-run (start, end) pairs per sample. Must be called before RLE
   multiply methods.

.. method:: ThreadingInstructions.right_multiply_rle(X)

   Compute G @ X via RLE. Accepts 1D array/list ``(num_sites,)`` or 2D array
   ``(num_sites, k)``. Returns numpy array of matching shape.

.. method:: ThreadingInstructions.left_multiply_rle(X)

   Compute G.T @ X via RLE. Accepts 1D array/list ``(num_samples,)`` or 2D array
   ``(num_samples, k)``. Returns numpy array of matching shape.

Dense Multiply
--------------

Falls back to materializing the full genotype matrix. O(n * m) per call.

.. method:: ThreadingInstructions.left_multiply(x, diploid=False, normalize=False)

   Compute G.T @ x (dense). If ``diploid``, uses diploid dosage matrix. If
   ``normalize``, uses mean-centered variance-normalized matrix.

.. method:: ThreadingInstructions.right_multiply(x, diploid=False, normalize=False)

   Compute G @ x (dense). Same options as ``left_multiply``.
