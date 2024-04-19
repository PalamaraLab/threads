import arg_needle_lib
import numpy as np
import logging
import time
import importlib
import os
from cyvcf2 import VCF

def mapping_string(carrier_sets, edges):
    if len(edges) == 0:
        return "NaN"
    elif len(edges) == 1:
        return f"-1,{edges[0].child.height:.4f},{edges[0].parent.height:.4f}"
    else:
        return ";".join([f"{'.'.join([str(c) for c in carrier_set])},{edge.child.height:.4f},{edge.parent.height:.4f}" for carrier_set, edge in zip(carrier_sets, edges)])

def get_leaves(arg, edge, position):
    leaves = []
    populate_leaves(arg, edge, position - arg.offset, leaves)
    return leaves

def populate_leaves(arg, edge, position, leaf_list):
    child = edge.child
    if arg.is_leaf(child.ID):
        return leaf_list.append(child.ID)
    else:
        # leaves = []
        for edge in child.child_edges_at(position):
            populate_leaves(arg, edge, position, leaf_list)

def map_region(argn, input, region, maf):
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

    # initialize counters etc
    maf_threshold = maf
    all_mappings = []
    n_attempted = 0
    n_mapped = 0
    n_parsimoniously_mapped = 0

    # iterate over VCF here
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
        _, mapping = arg_needle_lib.map_genotype_to_ARG_relate(arg, hap, float(pos - arg.offset), maf_threshold=maf)
        map_time += time.time() - mt

        if len(mapping) > 0:
            n_mapped += 1
        else:
            continue

        if len(mapping) == 1:
            all_mappings.append((name, pos, flipped, [[-1]],  mapping))
        else:
            all_mappings.append((name, pos, flipped, [get_leaves(arg, edge, pos) for edge in mapping],  mapping))

    n_mapped = sum(1 for m in all_mappings if len(m[4]) > 0)


    n_relate_mapped = n_mapped - n_parsimoniously_mapped
    return_strings = []
    for name, pos, flipped, carrier_sets, edges in all_mappings:
        return_strings.append(f"{name}\t{pos}\t{int(flipped)}\t{mapping_string(carrier_sets, edges)}\n")
    
    end_time = time.time()
    local_logger.info(f"Done region {region} in {end_time - start_time:.2f} (s)")
    return return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped
