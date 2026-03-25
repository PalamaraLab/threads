Output Formats
==============

VCF
---

.. module:: threads_arg.threads_to_vcf

.. autofunction:: threads_to_vcf

The C++ VCF writer (``VCFWriter``) uses buffered string output---each site is
assembled into a single string and flushed with one ``write()`` call.

Pgen
----

.. module:: threads_arg.threads_to_pgen

.. autofunction:: threads_to_pgen

Writes phased plink2 files (pgen/pvar/psam) using ``pgenlib.PgenWriter`` with
``all_phased=True``.
