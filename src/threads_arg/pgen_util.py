# import pgenlib
# import numpy as np

# def iterate_pgen(pgen, callback, start_idx=None, end_idx=None, **kwargs):
#     # Initialize read batching
#     num_samples = reader.get_raw_sample_ct()
#     num_sites = reader.get_variant_ct()
#     if start_idx is None:
#         start_idx = 0
#     if end_idx is None:
#         end_idx = num_sites
#     M = end_idx - start_idx

#     BATCH_SIZE = int(4e7 // num_samples)
#     n_batches = int(np.ceil(M / BATCH_SIZE))

#     logger.info("Finding singletons")
#     # Get singleton filter for the matching step
#     alleles_out = None
#     phased_out = None
#     i = 0
#     for b in range(n_batches):
#         b_start = b * BATCH_SIZE + start_idx
#         b_end = min(end_idx, (b+1) * BATCH_SIZE)
#         g_size = b_end - b_start
#         alleles_out = np.empty((g_size, 2 * num_samples), dtype=np.int32)
#         phased_out = np.empty((g_size, num_samples), dtype=np.uint8)
#         if np.any(phased_out == 0, axis=1):
#             raise RuntimeError("Unphased variants are currently not supported.")
#         reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
#         for g in alleles_out:
#             callback(i, g, **kwargs)
#             i += 1