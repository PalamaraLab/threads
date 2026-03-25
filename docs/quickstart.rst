Quick Start
===========

Installation
------------

.. code-block:: bash

   pip install threads-arg

For imputation support (requires pandas, scipy, cyvcf2):

.. code-block:: bash

   pip install 'threads-arg[impute]'

Inferring an ARG
----------------

.. code-block:: bash

   threads infer \
       --pgen input.pgen \
       --demography demography.txt \
       --mutation_rate 1.4e-8 \
       --mode wgs \
       --out output

This produces an ``output.threads`` file containing threading instructions.

Loading and using instructions
------------------------------

.. code-block:: python

   from threads_arg.serialization import load_instructions

   instructions = load_instructions("output.threads")

   # Basic properties
   print(instructions.num_samples)   # number of haploid samples
   print(instructions.num_sites)     # number of variant sites
   print(instructions.positions)     # site positions (bp)

   # Allele frequency spectrum (50-70x faster than arg-needle-lib)
   afs = instructions.allele_frequency_spectrum()

   # Total branch volume
   vol = instructions.total_volume()

   # Genotype matrix as numpy array
   geno = instructions.genotype_matrix_numpy()  # (n_sites, n_samples) int8

   # Generate mutations proportional to branch volume
   positions, genotypes = instructions.generate_mutations(mu=1e-8, seed=42)

   # Association testing
   pos, chi2, pval, ndesc = instructions.association_diploid(
       phenotypes, p_threshold=5e-8)

Genotype-matrix multiplication
------------------------------

Multiply without materializing the full genotype matrix:

.. code-block:: python

   import numpy as np

   # Tree shuttle (preferred for large m)
   instructions.prepare_tree_multiply()
   y = instructions.right_multiply_tree(x)   # G @ x
   z = instructions.left_multiply_tree(x)    # G.T @ x

   # Batch multiply (k vectors at once)
   Y = instructions.right_multiply_tree_batch(X_flat, k)

Export to standard formats
--------------------------

.. code-block:: bash

   # VCF to stdout
   threads vcf --threads output.threads

   # Plink2 pgen/pvar/psam
   threads pgen --threads output.threads --out prefix
