# This file is part of the Threads software suite.
# Copyright (C) 2024-2025 Threads Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import time
import logging
import arg_needle_lib
import numpy as np

logger = logging.getLogger(__name__)


def phase_distance(arg, arg_pos, target, het_carriers, hom_carriers):
    """
    Compute 'phasing distance' as described in the thesis.
    """
    if len(het_carriers) + 2 * len(hom_carriers) == 1:
        # negative so lower number (with longer pendant branch) wins
        return -arg.node(target).parent_edge_at(arg_pos).parent.height
    phase_distance = 0
    for carrier in het_carriers:
        phase_distance += min(arg.mrca(target, 2 * carrier, arg_pos).height, arg.mrca(target, 2 * carrier + 1, arg_pos).height)
    for hom in hom_carriers:
        phase_distance += arg.mrca(target, 2 * hom, arg_pos).height + arg.mrca(target, 2 * hom + 1, arg_pos).height
    return phase_distance


def threads_phase(scaffold, argn, ts, unphased, out):
    """
    Use an imputed arg to phase. Other input follows same shape as in SHAPEIT5-rare.
    """
    logger.info("Starting Threads-phase.")
    logger.info("WARNING: Threads-phase is experimental functionality.")
    start_time = time.time()
    unphased_vcf = VCF(unphased)
    scaffold_vcf = VCF(scaffold)
    true_vcf = VCF(unphased)
    phased_writer = Writer(out, true_vcf)
    if argn is None and ts is None:
        raise ValueError("Need either --argn or --ts")
    logger.info("Reading ARG...")
    try:
        arg = arg_needle_lib.deserialize_arg(argn)
    except:
        import tskit
        treeseq = tskit.load(ts)
        arg = arg_needle_lib.tskit_to_arg(treeseq)
    arg.populate_children_and_roots()

    i = 0
    logger.info("Phasing...")
    scaffold_empty = False
    v_scaffold = next(scaffold_vcf)
    s_scaffold = 0

    num_hets_found = 0
    # Main phasing routine
    for v in unphased_vcf:
        if not scaffold_empty and v.ID == v_scaffold.ID:
            # If variant exists in scaffold, just copy it
            v.genotypes = v_scaffold.genotypes
            s_scaffold += 1
            try:
                v_scaffold = next(scaffold_vcf)
            except StopIteration:
                scaffold_empty = True
        else:
            # Otherwise, do ARG-based phasing
            G = np.array(v.genotypes)
            g0, g1 = G[:, 0], G[:, 1]
            het_carriers = ((g0 + g1) == 1).nonzero()[0]
            hom_carriers = ((g0 + g1) == 2).nonzero()[0]
            arg_pos = max(v.end - arg.offset, 0)
            arg_pos = min(arg_pos, arg.end - 1)

            for target in het_carriers:
                num_hets_found += 1
                phase_distance_0 = phase_distance(arg, arg_pos, 2 * target, het_carriers, hom_carriers)
                phase_distance_1 = phase_distance(arg, arg_pos, 2 * target + 1, het_carriers, hom_carriers)
                if phase_distance_0 <= phase_distance_1:
                    G[target] = [1, 0, True]
                else:
                    G[target] = [0, 1, True]
            v.genotypes = G

        v.genotypes = v.genotypes
        phased_writer.write_record(v)
    logger.info(f"Done, in {time.time() - start_time} seconds")
    unphased_vcf.close()
    scaffold_vcf.close()
    phased_writer.close()
