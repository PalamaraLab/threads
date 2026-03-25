ThreadingInstructions
=====================

.. module:: threads_arg

The central data structure. Stores *n* haploid samples over *m* variant sites as
compact threading instructions: for each sample, a sequence of segments
(start_bp, target, tmrca) plus a list of mismatch sites.

.. autoclass:: ThreadingInstructions
   :members:
   :undoc-members:
   :show-inheritance:

Properties
----------

.. attribute:: ThreadingInstructions.num_samples

   Number of haploid samples (read-only int).

.. attribute:: ThreadingInstructions.num_sites

   Number of variant sites (read-only int).

.. attribute:: ThreadingInstructions.positions

   List of site positions in base pairs (read-only).

.. attribute:: ThreadingInstructions.start

   Region start in base pairs (read-only int).

.. attribute:: ThreadingInstructions.end

   Region end in base pairs (read-only int).

Accessors
---------

.. method:: ThreadingInstructions.all_starts()

   Return all segment start positions (bp) as a flat list, concatenated across samples.

.. method:: ThreadingInstructions.all_tmrcas()

   Return all segment TMRCAs (generations) as a flat list, concatenated across samples.

.. method:: ThreadingInstructions.all_targets()

   Return all segment target sample indices as a flat list, concatenated across samples.

.. method:: ThreadingInstructions.all_mismatches()

   Return all mismatch site indices as a flat list, concatenated across samples.

.. method:: ThreadingInstructions.sub_range(start_bp, end_bp)

   Return a new ThreadingInstructions restricted to sites in [start_bp, end_bp].

Genotype Matrix
---------------

.. method:: ThreadingInstructions.genotype_matrix_numpy()

   Return the full genotype matrix as a ``(num_sites, num_samples)`` int8 numpy
   array. Materializes and caches the matrix on first call.

.. method:: ThreadingInstructions.materialize_genotypes()

   Cache the dense genotype matrix in memory. Called automatically by multiply
   methods and ``genotype_matrix_numpy()``.

Variant Operations
------------------

.. method:: ThreadingInstructions.add_variants(new_positions, new_genotypes)

   Add variants to existing instructions without re-inference. ``new_genotypes``
   is a ``(n_new_sites, num_samples)`` int32 array.

.. method:: ThreadingInstructions.simplify_for_multiply()

   Strip TMRCAs and merge consecutive segments with the same target. Reduces
   multiply cost without changing genotypes.

.. method:: ThreadingInstructions.coarsen(min_sites)

   Merge segments shorter than ``min_sites`` (lossy). Reduces tree count for
   faster tree multiply at the cost of genotype accuracy.
