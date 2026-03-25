Command-Line Interface
======================

Threads is invoked as ``threads <command>``.

.. code-block:: text

   threads infer        Infer an ARG from genotype data
   threads vcf          Print genotypes to stdout in VCF format
   threads pgen         Write genotypes to pgen/pvar/psam
   threads convert      Convert threading instructions to arg-needle format
   threads allele_ages  Estimate allele ages
   threads map          Map mutations to ARG branches
   threads impute       Impute ungenotyped variants

infer
-----

.. code-block:: bash

   threads infer \
       --pgen input.pgen \
       --demography demography.txt \
       --mutation_rate 1.4e-8 \
       --mode {array,wgs} \
       --out output \
       [--map genetic_map.txt] \
       [--recombination_rate 1.3e-8] \
       [--fit_to_data] \
       [--normalize] \
       [--num_threads 4] \
       [--save_metadata] \
       [--region chr1:1000000-2000000]

vcf
---

.. code-block:: bash

   threads vcf --threads output.threads [--variants input.pvar] [--samples samples.txt]

Writes phased VCF (``0|0``, ``0|1``, etc.) to stdout.

pgen
----

.. code-block:: bash

   threads pgen --threads output.threads --out prefix [--variants input.pvar] [--samples samples.txt]

Writes ``prefix.pgen``, ``prefix.pvar``, ``prefix.psam``.
