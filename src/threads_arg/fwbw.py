# Code adapted from an implementation of the Li-Stephens algorithm
# available at: https://github.com/astheeggeggs/lshmm
import numpy as np
import logging

from threads_arg import forwards_ls_hap, backwards_ls_hap

logger = logging.getLogger(__name__)

MISSING = -9


def checks(reference_panel, query, mutation_rate, recombination_rates):
    ref_shape = reference_panel.shape
    ploidy = len(ref_shape) - 1
    assert ploidy == 1
    m, _ = ref_shape

    if not (query.shape[1] == m):
        raise ValueError(
            "Number of variants in query does not match reference_panel. Please ensure variant x sample matrices are passed."
        )

    # Ensure that the mutation rate is a scalar
    if not isinstance(mutation_rate, float):
        raise ValueError("mutation_rate is not a float")

    # Ensure that the recombination probabilities is a vector of length m
    if recombination_rates.shape[0] != m:
        raise ValueError(f"recombination_rates are not a vector of length m: {m}")

    return


def set_emission_probabilities(reference_panel, query, mutation_rate):
    m, n = reference_panel.shape

    # Vectorized check: a site is biallelic if it has both 0s and 1s across
    # the combined reference+query panel. For binary data this is equivalent
    # to checking that the row min != row max (after including the query).
    combined_min = np.minimum(reference_panel.min(axis=1), query.min(axis=0))
    combined_max = np.maximum(reference_panel.max(axis=1), query.max(axis=0))
    is_polymorphic = combined_min != combined_max
    n_alleles = np.where(is_polymorphic, np.int8(2), np.int8(1))

    if not np.all((n_alleles == 2) | (n_alleles == 1)):
        raise ValueError("Only fixed or bi-allelic sites allowed")

    if mutation_rate is None:
        # Set the mutation rate to be the proposed mutation rate in Li and Stephens (2003).
        theta_tilde = 1 / np.sum([1 / k for k in range(1, n - 1)])
        mutation_rate = 0.5 * (theta_tilde / (n + theta_tilde))

    mutation_rate = mutation_rate * np.ones(m)

    # Evaluate emission probabilities here, using the mutation rate
    e = np.zeros((m, 2))
    e[:, 0] = mutation_rate
    e[:, 1] = 1 - mutation_rate
    # Invariant sites: emission for mismatch is 0, match is 1
    invariant = n_alleles == 1
    e[invariant, 0] = 0
    e[invariant, 1] = 1
    return e


def fwbw(reference_panel,
         query,
         recombination_rates,
         mutation_rate):
    # Input checks
    checks(reference_panel, query, mutation_rate, recombination_rates)
    m, n = reference_panel.shape

    # Get emissions
    emissions = set_emission_probabilities(reference_panel, query, mutation_rate)

    # Run forwards (C++)
    forward_array, fwd_norm_factor = forwards_ls_hap(
        reference_panel.astype(np.float64),
        query.ravel().astype(np.float64),
        emissions,
        recombination_rates)

    # Run backwards (C++)
    backward_array = backwards_ls_hap(
        reference_panel.astype(np.float64),
        query.ravel().astype(np.float64),
        emissions,
        fwd_norm_factor,
        recombination_rates)

    # Return posterior
    return forward_array * backward_array
