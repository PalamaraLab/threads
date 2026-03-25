Genotype-Matrix Multiplication
==============================

Multiply the implicit genotype matrix by vectors without materializing the full
matrix. Two strategies are available: tree shuttle (preferred) and RLE.

Tree Shuttle
------------

Propagates values through the interval-tree structure. Preferred for large *m*
because it replaces the O(m) sweep with O(n * n_intervals), where n_intervals is
typically 10-100x smaller than m.

**Complexity:** O(n * n_intervals + mismatches) per multiply.

.. method:: ThreadingInstructions.prepare_tree_multiply()

   Precompute interval-tree structure. Must be called before any tree multiply
   method. Uses reference-counted genotype cache with O(m * tree_depth) peak
   memory.

.. method:: ThreadingInstructions.right_multiply_tree(x)

   Compute G @ x. x has length ``num_sites``. Returns list of length
   ``num_samples``.

.. method:: ThreadingInstructions.left_multiply_tree(x)

   Compute G.T @ x. x has length ``num_samples``. Returns list of length
   ``num_sites``.

.. method:: ThreadingInstructions.right_multiply_tree_batch(x_flat, k)

   Compute G @ X for k vectors. ``x_flat`` has length ``num_sites * k``.
   Returns flat list of length ``num_samples * k``.

.. method:: ThreadingInstructions.left_multiply_tree_batch(x_flat, k)

   Compute G.T @ X for k vectors. ``x_flat`` has length ``num_samples * k``.
   Returns flat list of length ``num_sites * k``.

.. method:: ThreadingInstructions.right_multiply_tree_batch_numpy(x_flat, k)

   Same as ``right_multiply_tree_batch`` but accepts and returns numpy arrays.

.. method:: ThreadingInstructions.left_multiply_tree_batch_numpy(x_flat, k)

   Same as ``left_multiply_tree_batch`` but accepts and returns numpy arrays.

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

.. method:: ThreadingInstructions.right_multiply_rle(x)

   Compute G @ x via RLE. x has length ``num_sites``.

.. method:: ThreadingInstructions.left_multiply_rle(x)

   Compute G.T @ x via RLE. x has length ``num_samples``.

.. method:: ThreadingInstructions.right_multiply_rle_batch(x_flat, k)

   Compute G @ X via RLE for k vectors.

.. method:: ThreadingInstructions.left_multiply_rle_batch(x_flat, k)

   Compute G.T @ X via RLE for k vectors.

Dense Multiply
--------------

Falls back to materializing the full genotype matrix. O(n * m) per call.

.. method:: ThreadingInstructions.left_multiply(x, diploid=False, normalize=False)

   Compute G.T @ x (dense). If ``diploid``, uses diploid dosage matrix. If
   ``normalize``, uses mean-centered variance-normalized matrix.

.. method:: ThreadingInstructions.right_multiply(x, diploid=False, normalize=False)

   Compute G @ x (dense). Same options as ``left_multiply``.
