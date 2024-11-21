# import threads_arg
# import pgenlib
# import numpy as np

# def estimate_allele_ages(instructions, physical_positions, pgen, start_idx, end_idx):
#     # NB! this code is largely copied from the infer() function, 
#     # consider writing a batched readier that does this consistently
#     age_estimator = threads_arg.AgeEstimator(instructions)

#     reader = pgenlib.PgenReader(pgen.encode())
#     num_samples = reader.get_raw_sample_ct()
#     num_sites = reader.get_variant_ct()
#     assert physical_positions == num_sites
#     site_indicator = (start <= physical_positions) & (physical_positions <= end)

#     # Initialize read batching
#     M = len(physical_positions)
#     assert num_sites == M
#     BATCH_SIZE = int(4e7 // num_samples)
#     n_batches = int(np.ceil(M / BATCH_SIZE))

#     logger.info("Finding singletons")
#     # Get singleton filter for the matching step
#     alleles_out = None
#     phased_out = None
#     site_counter = 0
#     for b in range(n_batches):
#         b_start = b * BATCH_SIZE
#         b_end = min(M, (b+1) * BATCH_SIZE)
#         g_size = b_end - b_start
#         alleles_out = np.empty((g_size, 2 * num_samples), dtype=np.int32)
#         phased_out = np.empty((g_size, num_samples), dtype=np.uint8)
#         reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
#         # Filter out unphased variants and singletons
#         for g in alleles_out:
#             age_estimator.process_site(g)
#     return age_estimator.get_inferred_ages
