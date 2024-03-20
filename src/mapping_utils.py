# import ray
# import click
import arg_needle_lib
import numpy as np
import logging
import time
import importlib
import os
from cyvcf2 import VCF

# os.environ["RAY_DEDUP_LOGS"] = "0"

# logging.basicConfig(
#     format='%(asctime)s %(levelname)-8s %(message)s',
#     level=logging.INFO,
#     datefmt='%Y-%m-%d %H:%M:%S')

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

# # see https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
# def split(a, n):
#     k, m = divmod(len(a), n)
#     return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

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
    # ret_strings = []
    return_strings = []
    for name, pos, flipped, carrier_sets, edges in all_mappings:
        return_strings.append(f"{name}\t{pos}\t{int(flipped)}\t{mapping_string(carrier_sets, edges)}\n")
    
    end_time = time.time()
    local_logger.info(f"Done region {region} in {end_time - start_time:.2f} (s)")
    return return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped


# @click.command()
# @click.option("--argn", help="Path to .argn file with results")
# @click.option("--out", help="Where to save results")
# @click.option("--maf", type=float, default=0.02, help="Don't map stuff with MAF above this")
# @click.option("--input", type=str, help="Path to bcf/vcf with genotypes to map. Most have AC/AN fields")
# @click.option("--region", type=str, help="Of format chr:start-end (both inclusive)")
# @click.option("--threads", type=int, help="Of format chr:start-end (both inclusive)", default=1)
# def map_mutations_to_arg(argn, out, maf, input, region, threads):
#     logging.info("Starting mutation-mapping with the following arguments:")
#     logging.info(f"argn:   {argn}")
#     logging.info(f"out:    {out}")
#     logging.info(f"maf:    {maf}")
#     logging.info(f"input:  {input}")
#     logging.info(f"region: {region}")
#     logging.info(f"threads: {threads}")
#     start_time = time.time()

#     actual_num_threads = min(len(os.sched_getaffinity(0)), threads)
#     logging.info(f"Requested {threads} threads, found {actual_num_threads}.")

#     return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped = None, None, None, None
#     if actual_num_threads == 1:
#         return_strings, n_attempted, n_parsimoniously_mapped, n_relate_mapped = map_region(argn, input, region, maf)
#     else:
#         logging.info("Parsing VCF")
#         vcf = VCF(input)
#         positions = [record.POS for record in vcf(region)]
#         assert len(vcf.seqnames) == 1
#         contig = vcf.seqnames[0]

#         # split into subregions
#         split_positions = split(positions, actual_num_threads)
#         subregions = [f"{contig}:{pos[0]}-{pos[-1]}" for pos in split_positions]
#         ray.init()
#         map_region_remote = ray.remote(map_region)
#         results = ray.get([map_region_remote.remote(
#             argn, input, subregion, maf
#         ) for subregion in subregions])
#         ray.shutdown()
#         return_strings = []
#         n_attempted, n_parsimoniously_mapped, n_relate_mapped = 0, 0, 0
#         for rets, natt, npars, nrel in results:
#             return_strings += rets
#             n_attempted += natt
#             n_parsimoniously_mapped += npars
#             n_relate_mapped += nrel

#     end_time = time.time()
#     logging.info(f"Done in (s): {end_time-start_time:.3f}")
#     n_mapped = n_relate_mapped + n_parsimoniously_mapped
#     n_unsuccessful = n_attempted - n_mapped

#     logging.info(f"Attempted to map {n_attempted} variants.")
#     logging.info(f"{n_parsimoniously_mapped} ({100 * n_parsimoniously_mapped / n_attempted:.3f}%) had MAC≤4 and mapped trivially.")
#     logging.info(f"{n_relate_mapped} ({100 * n_relate_mapped / n_attempted:.3f}%) MAC>4 variants mapped.")
#     logging.info(f"{n_unsuccessful} ({100 * n_unsuccessful / n_attempted:.3f}%) variants did not map.")
#     logging.info(f"Writing mutation mappings to {out}")

#     # Output has columns 
#     # variant_id: string
#     # pos: int (base-pairs)
#     # flipped: bool
#     # mapping_string of the following format:
#     #   If uniquely mapped, a string "-1,edge_lower,edge_upper"
#     #   If multiply mapped, strings "leaf_1...leaf_k1,edge1_lower,edge1_upper;...;leaf_1...leafkN,edgeN_lower,edgeN_upper"
#     #       telling us the range (in generations) of each edge and the leaves it subtends (-1 if all carriers)
#     with open(out, "w") as outfile:
#         for string in return_strings:
#             outfile.write(string)
#     logging.info(f"Total runtime {time.time() - start_time:.2f}")

# if __name__ == "__main__":
#     map_mutations_to_arg()
