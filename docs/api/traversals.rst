ARG Traversals and Statistics
=============================

These methods reconstruct the implicit coalescent tree from threading
instructions on the fly using union-find, without materializing an explicit ARG.

Tree Traversals
---------------

.. method:: ThreadingInstructions.visit_clades(callback)

   Traverse all clades in the implicit ARG. For each clade, calls
   ``callback(descendants, height, start_bp, end_bp)`` where ``descendants``
   is a list of haploid sample indices.

   Uses the full union-find with member tracking.

.. method:: ThreadingInstructions.visit_branches(callback)

   Traverse all branches in the implicit ARG. For each branch, calls
   ``callback(descendants, child_height, parent_height, start_bp, end_bp)``.

   Uses the full union-find with member tracking.

Summary Statistics
------------------

.. method:: ThreadingInstructions.total_volume()

   Total branch volume: sum of (parent_height - child_height) * span over all
   branches. Uses a lightweight union-find tracking only heights (no member
   lists).

   **Complexity:** O(n * L * alpha(n)) where L is the number of distinct tree
   intervals.

   :returns: Total volume as a float.

.. method:: ThreadingInstructions.allele_frequency_spectrum()

   Branch-length allele frequency spectrum: ``afs[k]`` is the total volume of
   branches with exactly *k* descendants. Uses a count-only union-find (sizes
   but no member lists), ~50-70x faster than arg-needle-lib's
   ``bitset_volume_map``.

   **Complexity:** O(n * L * alpha(n)), same as total_volume.

   :returns: List of floats, length ``num_samples + 1``.

Association Testing
-------------------

.. method:: ThreadingInstructions.association_diploid(phenotypes, p_threshold=5e-8)

   Clade-based association testing. For each clade, computes diploid dosage and
   tests against phenotypes via chi-square (Pearson r^2 * n_dip, 1 df).

   :param phenotypes: List of phenotype values, length ``num_samples / 2``.
   :param p_threshold: Only report clades with p-value below this threshold.
   :returns: Tuple of ``(positions, chi2_stats, p_values, n_descendants)``.

Mutation Generation
-------------------

.. method:: ThreadingInstructions.generate_mutations(mutation_rate, seed=42)

   Poisson-sample mutations on branches proportional to volume. Each mutation is
   placed uniformly within its branch's genomic span and inherits the branch's
   descendant set as carrier genotypes.

   :param mutation_rate: Per-base per-generation mutation rate.
   :param seed: Random seed.
   :returns: Tuple of ``(positions, genotypes)`` where positions is an int32
      array of length n_mutations and genotypes is an int32 array of shape
      ``(n_mutations, num_samples)``.
