# Code adapted from an implementation of the Li-Stephens algorithm
# available at: https://github.com/astheeggeggs/lshmm
import numba
import numpy as np
import logging
import os

logger = logging.getLogger(__name__)

_DISABLE_NUMBA = os.environ.get("LSHMM_DISABLE_NUMBA", "0")

try:
    ENABLE_NUMBA = {"0": True, "1": False}[_DISABLE_NUMBA]
except KeyError as e:
    raise KeyError(
        "Environment variable 'LSHMM_DISABLE_NUMBA' must be '0' or '1'"
    ) from e

if not ENABLE_NUMBA:
    logger.warning(
        "Numba globally disabled, performance will be drastically reduced."
    )


DEFAULT_NUMBA_ARGS = {
    "nopython": True,
    "cache": True,
}


def numba_njit(func, **kwargs):
    if ENABLE_NUMBA:
        return numba.jit(func, **{**DEFAULT_NUMBA_ARGS, **kwargs})
    else:
        return func

MISSING = -9

@numba_njit
def forwards_ls_hap(n, m, H, s, e, r, norm=True):
    """Matrix based haploid LS forward algorithm using numpy vectorisation."""
    # Initialise
    F = np.zeros((m, n))
    r_n = r / n

    if norm:
        c = np.zeros(m)
        for i in range(n):
            F[0, i] = (
                1 / n * e[0, np.int64(np.equal(H[0, i], s[0, 0]) or s[0, 0] == MISSING)]
            )
            c[0] += F[0, i]

        for i in range(n):
            F[0, i] *= 1 / c[0]

        # Forwards pass
        for l in range(1, m):
            for i in range(n):
                F[l, i] = F[l - 1, i] * (1 - r[l]) + r_n[l]
                F[l, i] *= e[
                    l, np.int64(np.equal(H[l, i], s[0, l]) or s[0, l] == MISSING)
                ]
                c[l] += F[l, i]

            for i in range(n):
                F[l, i] *= 1 / c[l]
        # Log-likelihood: ll = np.sum(np.log10(c))
    else:
        c = np.ones(m)
        for i in range(n):
            F[0, i] = (
                1 / n * e[0, np.int64(np.equal(H[0, i], s[0, 0]) or s[0, 0] == MISSING)]
            )

        # Forwards pass
        for l in range(1, m):
            for i in range(n):
                F[l, i] = F[l - 1, i] * (1 - r[l]) + np.sum(F[l - 1, :]) * r_n[l]
                F[l, i] *= e[
                    l, np.int64(np.equal(H[l, i], s[0, l]) or s[0, l] == MISSING)
                ]
        # Log-likelihood: ll = np.log10(np.sum(F[m - 1, :]))
    return F, c

@numba_njit
def backwards_ls_hap(n, m, H, s, e, c, r):
    """Matrix based haploid LS backward algorithm using numpy vectorisation."""
    # Initialise
    B = np.zeros((m, n))
    for i in range(n):
        B[m - 1, i] = 1
    r_n = r / n

    # Backwards pass
    for l in range(m - 2, -1, -1):
        tmp_B = np.zeros(n)
        tmp_B_sum = 0
        for i in range(n):
            tmp_B[i] = (
                e[
                    l + 1,
                    np.int64(
                        np.equal(H[l + 1, i], s[0, l + 1]) or s[0, l + 1] == MISSING
                    ),
                ]
                * B[l + 1, i]
            )
            tmp_B_sum += tmp_B[i]
        for i in range(n):
            B[l, i] = r_n[l + 1] * tmp_B_sum
            B[l, i] += (1 - r[l + 1]) * tmp_B[i]
            B[l, i] *= 1 / c[l + 1]

    return B


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
    n_alleles = np.int8(
        [
            len(np.unique(np.append(reference_panel[j, :], query[:, j])))
            for j in range(reference_panel.shape[0])
        ]
    )

    if not np.all((n_alleles == 2) | (n_alleles == 1)):
        raise ValueError("Only fixed or bi-allelic sites allowed")

    if mutation_rate is None:
        # Set the mutation rate to be the proposed mutation rate in Li and Stephens (2003).
        theta_tilde = 1 / np.sum([1 / k for k in range(1, n - 1)])
        mutation_rate = 0.5 * (theta_tilde / (n + theta_tilde))

    mutation_rate = mutation_rate * np.ones(m)

    # Evaluate emission probabilities here, using the mutation rate
    e = np.zeros((m, 2))
    for j in range(m):
        if n_alleles[j] == 1:  # In case we're at an invariant site
            e[j, 0] = 0
            e[j, 1] = 1
        else:
            e[j, 0] = mutation_rate[j]
            e[j, 1] = 1 - mutation_rate[j]
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

    # Run forwards
    forward_array, fwd_norm_factor = forwards_ls_hap(
        n, m, reference_panel, query, emissions, recombination_rates, norm=True)

    # Run backwards
    backward_array = backwards_ls_hap(
        n, m, reference_panel, query, emissions, fwd_norm_factor, recombination_rates)

    # Return posterior
    return forward_array * backward_array
