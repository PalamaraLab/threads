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

import argparse
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)


def goodbye():
    # Credit to a nameless contribution at https://www.asciiart.eu/miscellaneous/dna
    print(
    """
    `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
       `=`,'=/     `-` '-'     `-` '-'     `=`,'=/
         y==/  Thank you for using Threads!  y==/
       ,=,-<=`.    ,=,-,=-.    ,=,-,=-.    ,=,-<=`.
    ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_""")


def main():
    parser = argparse.ArgumentParser(prog="threads")
    subparsers = parser.add_subparsers(dest="command")

    # infer
    p_infer = subparsers.add_parser("infer", help="Infer an ARG from genotype data")
    p_infer.add_argument("--pgen", required=True, help="Path to input genotypes in pgen format")
    p_infer.add_argument("--map", help="Path to genotype map in SHAPEIT format")
    p_infer.add_argument("--recombination_rate", default=1.3e-8, type=float, help="Genome-wide recombination rate. Ignored if a map is passed")
    p_infer.add_argument("--demography", required=True, help="Path to input genotype")
    p_infer.add_argument("--mode", required=True, choices=["array", "wgs"], default="wgs", help="Inference mode (wgs or array)")
    p_infer.add_argument("--fit_to_data", action="store_true", default=False, help="If specified, Threads performs a post-processing step to ensure the inferred ARG contains an edge matching each input mutation.")
    p_infer.add_argument("--normalize", action="store_true", default=False, help="If specified, Threads will normalize output to fit with the input demography.")
    p_infer.add_argument("--allele_ages", default=None, help="Allele ages used for post-processing with the --fit_to_data flag, otherwise ignored. If not specified, allele ages are inferred automatically.")
    p_infer.add_argument("--query_interval", type=float, default=0.01, help="Hyperparameter for the preliminary haplotype matching in cM")
    p_infer.add_argument("--match_group_interval", type=float, default=0.5, help="Hyperparameter for the preliminary haplotype matching in cM")
    p_infer.add_argument("--mutation_rate", required=True, type=float, default=1.4e-8, help="Genome-wide mutation rate")
    p_infer.add_argument("--num_threads", type=int, default=1, help="Number of computational threads to request")
    p_infer.add_argument("--region", help="Region of genome in chr:start-end format for which ARG is output. The full genotype is still used for inference")
    p_infer.add_argument("--max_sample_batch_size", help="Max number of LS processes run simultaneously per thread", default=None, type=int)
    p_infer.add_argument("--save_metadata", action="store_true", default=False, help="If specified, the output will include sample/variant metadata (sample IDs, marker names, allele symbols, etc).")
    p_infer.add_argument("--out")

    # convert
    p_convert = subparsers.add_parser("convert", help="Convert Threads ARGs to ARG-Needle or tskit format")
    p_convert.add_argument("--threads", required=True, help="Path to an input .threads file")
    p_convert.add_argument("--argn", default=None, help="Path to an output .argn file")
    p_convert.add_argument("--tsz", default=None, help="Path to an output .tsz file")
    p_convert.add_argument("--add_mutations", action="store_true", default=False, help="If passed, mutations are parsimoniously added to the output ARG. This may result in a high number of mutations if the --fit_to_data flag was not used.")

    # allele_ages
    p_allele_ages = subparsers.add_parser("allele_ages", help="Infer allele ages from a Threads ARG")
    p_allele_ages.add_argument("--threads", required=True, help="Path to an input .threads file.")
    p_allele_ages.add_argument("--out", required=True, help="Path to output.")
    p_allele_ages.add_argument("--num_threads", type=int, help="Size of processor pool to process batches", default=None)

    # map
    p_map = subparsers.add_parser("map", help="Map genotypes to an ARG in ARG-Needle format")
    p_map.add_argument("--argn", help="Path to input .argn file")
    p_map.add_argument("--out", help="Path to output .mut file")
    p_map.add_argument("--maf", type=float, default=0.02, help="Do not store entries with MAF above this")
    p_map.add_argument("--input", type=str, help="Path to bcf/vcf with genotypes to map with AC/AN fields")
    p_map.add_argument("--region", type=str, help="Region in chr:start-end format (start and end inclusive)")
    p_map.add_argument("--num_threads", type=int, help="Number of computational threads to request", default=1)

    # impute
    p_impute = subparsers.add_parser("impute", help="Impute missing genotypes using a reference panel")
    p_impute.add_argument("--panel", required=True, help="pgen array panel")
    p_impute.add_argument("--target", required=True, help="pgen array targets")
    p_impute.add_argument("--mut", required=True, help="pgen array targets")
    p_impute.add_argument("--map", required=True, help="Path to genotype map in SHAPEIT format")
    p_impute.add_argument("--mutation_rate", type=float, help="Per-site-per-generation SNP mutation rate", default=1.4e-8)
    p_impute.add_argument("--demography", required=True, help="Path to file containing demographic history")
    p_impute.add_argument("--out", help="Path to output .vcf file", default=None)
    p_impute.add_argument("--stdout", help="Redirect output to stdout (will disable logging)", action="store_true")
    p_impute.add_argument("--region", required=True, type=str, help="Region in chr:start-end format (start and end inclusive)")

    # vcf
    p_vcf = subparsers.add_parser("vcf", help="Print genotypes from Threads ARGs to stdout in VCF format")
    p_vcf.add_argument("--threads", required=True, help="Path to input .threads file")
    p_vcf.add_argument("--variants", default=None, help="Path to .pvar or .bim file with variant information")
    p_vcf.add_argument("--samples", default=None, help="Path to a file with one sample ID per line")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        parser.exit(2)

    if args.command == "infer":
        from .infer import threads_infer
        kwargs = vars(args)
        del kwargs["command"]
        threads_infer(**kwargs)
        goodbye()

    elif args.command == "convert":
        from .convert import threads_convert
        kwargs = vars(args)
        del kwargs["command"]
        threads_convert(**kwargs)
        goodbye()

    elif args.command == "allele_ages":
        from .allele_ages import estimate_allele_ages
        kwargs = vars(args)
        del kwargs["command"]
        estimate_allele_ages(**kwargs)
        goodbye()

    elif args.command == "map":
        from .map_mutations_to_arg import threads_map_mutations_to_arg
        kwargs = vars(args)
        del kwargs["command"]
        threads_map_mutations_to_arg(**kwargs)
        goodbye()

    elif args.command == "impute":
        import sys
        if (args.stdout and args.out) or not (args.stdout or args.out):
            print("Either --out or --stdout must be specified", file=sys.stderr)
            exit(1)
        from .impute import Impute
        Impute(args.panel, args.target, args.map, args.mut, args.demography, args.out, args.region, args.mutation_rate)
        if not args.stdout:
            goodbye()

    elif args.command == "vcf":
        from .threads_to_vcf import threads_to_vcf
        threads_to_vcf(args.threads, samples=args.samples, variants=args.variants)


if __name__ == "__main__":
    main()
