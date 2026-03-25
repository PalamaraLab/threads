Utilities
=========

.. module:: threads_arg.utils

Pgen I/O
--------

.. autofunction:: read_all_genotypes

.. autofunction:: pgen_chunk_iterator

.. autofunction:: iterate_pgen

.. autofunction:: read_positions_and_ids

.. autofunction:: read_variant_metadata

.. autofunction:: read_sample_names

Genetic Maps and Demography
----------------------------

.. autofunction:: read_map_file

.. autofunction:: make_recombination_from_map_and_pgen

.. autofunction:: make_constant_recombination_from_pgen

.. autofunction:: parse_demography

.. autoclass:: VariantMetadata
   :members:
   :undoc-members:

Helpers
-------

.. autofunction:: split_list

.. autofunction:: parse_region_string

.. autofunction:: default_process_count

.. autoclass:: TimerTotal
   :members:
