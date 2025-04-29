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

import importlib
import os
import time
import logging

os.environ["RAY_DEDUP_LOGS"] = "0"
import ray
import numpy as np
import arg_needle_lib

from cyvcf2 import VCF
from .utils import split_list, default_process_count

logger = logging.getLogger(__name__)


def _mapping_string(carrier_sets, edges):
    if len(edges) == 0:
        return "NaN"
    elif len(edges) == 1:
        return f"-1,{edges[0].child.height:.4f},{edges[0].parent.height:.4f}"
    else:
        return ";".join([f"{'.'.join([str(c) for c in carrier_set])},{edge.child.height:.4f},{edge.parent.height:.4f}" for carrier_set, edge in zip(carrier_sets, edges)])


def get_leaf_ids_at(arg, edge, position):
    """
    Recurse into all child edges at position to get leave node IDs
    """
    leaves = []
    _recursive_get_leaf_ids_at(arg, edge, position - arg.offset, leaves)
    return leaves


def _recursive_get_leaf_ids_at(arg, edge, position, leaf_list):
    child = edge.child
    if arg.is_leaf(child.ID):
        return leaf_list.append(child.ID)
    else:
        for edge in child.child_edges_at(position):
            _recursive_get_leaf_ids_at(arg, edge, position, leaf_list)


def _map_region(argn, input, region, maf_threshold):
    logging.shutdown()
    importlib.reload(logging)
    pid = os.getpid()
    logging.basicConfig(format=f"%(asctime)s %(levelname)-8s PID {pid} %(message)s",
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')

    local_logger = logging.getLogger(__name__)
    start_time = time.time()
    local_logger.info(f"Starting region {region}...")
    arg = arg_needle_lib.deserialize_arg(argn)
    arg.populate_children_and_roots()

    # Initialize counters
    all_mappings = []
    n_attempted = 0
    n_mapped = 0
    n_parsimoniously_mapped = 0

    # Iterate over VCF records
    read_time = 0
    map_time = 0
    vcf = VCF(input)
    for record in vcf(region):
        ac = int(record.INFO.get("AC"))
        an = int(record.INFO.get("AN"))
        af = ac / an
        mac = min(ac, an - ac)
        maf = min(af, 1 - af)
        flipped = af > 0.5

        # Apply MAF filter
        if maf > maf_threshold or maf == 0:
            continue

        n_attempted += 1
        if mac <= 4:
            n_parsimoniously_mapped += 1

        # Thresholds passed, so we fetch genotype and attempt to map
        name = record.ID
        pos = record.POS

        rt = time.time()
        hap = np.array(record.genotypes)[:, :2].flatten()
        read_time += time.time() - rt
        assert len(hap) == len(arg.leaf_ids)
        if flipped:
            hap = 1 - hap

        mt = time.time()
        mapping, _ = arg_needle_lib.map_genotype_to_ARG_approximate(arg, hap, float(pos - arg.offset))
        map_time += time.time() - mt

        if len(mapping) > 0:
            n_mapped += 1
        else:
            continue

        if len(mapping) == 1:
            all_mappings.append((name, pos, flipped, [[-1]],  mapping))
        else:
            leaf_ids = [get_leaf_ids_at(arg, edge, pos) for edge in mapping]
            all_mappings.append((name, pos, flipped, leaf_ids, mapping))

    n_mapped = sum(1 for m in all_mappings if len(m[4]) > 0)


    n_relate_mapped = n_mapped - n_parsimoniously_mapped
    return_strings = []
    for name, pos, flipped, carrier_sets, edges in all_mappings:
        return_strings.append(f"{name}\t{pos}\t{int(flipped)}\t{_mapping_string(carrier_sets, edges)}\n")

    end_time = time.time()
    local_logger.info(f"Done region {region} in {end_time - start_time:.2f} (s)")
    return return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped


def threads_map_mutations_to_arg(argn, out, maf, input, region, num_threads):
    """
    Map mutations to an ARG using a method based on Speidel et al. (2019) and save output to a .mut file to inform imputation.

    Output has the following columns
    variant_id: string
    pos: int (base-pairs)
    flipped: bool
    mapping_string of the following format:
      If uniquely mapped, a string "-1,edge_lower,edge_upper"
      If multiply mapped, strings "leaf_1...leaf_k1,edge1_lower,edge1_upper;...;leaf_1...leafkN,edgeN_lower,edgeN_upper"
          telling us the range (in generations) of each edge and the leaves it subtends (-1 if all carriers)
    """
    logger.info(f"Starting Threads-map with parameters")
    logger.info(f"  argn:    {argn}")
    logger.info(f"  out:     {out}")
    logger.info(f"  maf:     {maf}")
    logger.info(f"  input:   {input}")
    logger.info(f"  region:  {region}")
    logger.info(f"  threads: {num_threads}")
    start_time = time.time()

    actual_num_threads = min(default_process_count(), num_threads)
    logger.info(f"Requested {num_threads} threads, found {actual_num_threads}.")

    return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped = None, None, None, None
    if actual_num_threads == 1:
        return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped = _map_region(argn, input, region, maf)
    else:
        logger.info("Parsing VCF")
        vcf = VCF(input)
        positions = [record.POS for record in vcf(region)]
        assert len(vcf.seqnames) == 1
        contig = vcf.seqnames[0]

        # split into subregions
        split_positions = split_list(positions, actual_num_threads)
        subregions = [f"{contig}:{pos[0]}-{pos[-1]}" for pos in split_positions]
        ray.init()
        map_region_remote = ray.remote(_map_region)
        results = ray.get([map_region_remote.remote(
            argn, input, subregion, maf
        ) for subregion in subregions])
        ray.shutdown()
        return_strings = []
        n_attempted, n_parsimoniously_mapped, n_relate_mapped = 0, 0, 0
        for rets, natt, npars, nrel in results:
            return_strings += rets
            n_attempted += natt
            n_parsimoniously_mapped += npars
            n_relate_mapped += nrel

    n_mapped = n_relate_mapped + n_parsimoniously_mapped
    n_unsuccessful = n_attempted - n_mapped

    logger.info(f"Attempted to map {n_attempted} variants.")
    logger.info(f"{n_parsimoniously_mapped} ({100 * n_parsimoniously_mapped / n_attempted:.3f}%) had MAC≤4 and mapped trivially.")
    logger.info(f"{n_relate_mapped} ({100 * n_relate_mapped / n_attempted:.3f}%) MAC>4 variants mapped.")
    logger.info(f"{n_unsuccessful} ({100 * n_unsuccessful / n_attempted:.3f}%) variants did not map.")
    logger.info(f"Writing mutation mappings to {out}")
    logger.info(f"Done in (s): {time.time()-start_time:.3f}")

    with open(out, "w") as outfile:
        for string in return_strings:
            outfile.write(string)
    logger.info(f"Total runtime {time.time() - start_time:.2f}")
