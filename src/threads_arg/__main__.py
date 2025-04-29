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

import click
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


@click.group()
def main():
    pass

@main.command()
@click.option("--pgen", required=True, help="Path to input genotypes in pgen format")
@click.option("--map", required=True, help="Path to genotype map in SHAPEIT format")
@click.option("--recombination_rate", default=1.3e-8, type=float, help="Genome-wide recombination rate. Ignored if a map is passed")
@click.option("--demography", required=True, help="Path to input genotype")
@click.option("--mode", required=True, type=click.Choice(['array', 'wgs']), default="wgs", help="Inference mode (wgs or array)")
@click.option("--fit_to_data", is_flag=True, default=False, help="If specified, Threads performs a post-processing step to ensure the inferred ARG contains an edge matching each input mutation.")
@click.option("--allele_ages", default=None, help="Allele ages used for post-processing with the --fit_to_data flag, otherwise ignored. If not specified, allele ages are inferred automatically.")
@click.option("--query_interval", type=float, default=0.01, help="Hyperparameter for the preliminary haplotype matching in cM")
@click.option("--match_group_interval", type=float, default=0.5, help="Hyperparameter for the preliminary haplotype matching in cM")
@click.option("--mutation_rate", required=True, type=float, default=1.4e-8, help="Genome-wide mutation rate")
@click.option("--num_threads", type=int, default=1, help="Number of computational threads to request")
@click.option("--region", help="Region of genome in chr:start-end format for which ARG is output. The full genotype is still used for inference")
@click.option("--max_sample_batch_size", help="Max number of LS processes run simultaneously per thread", default=None, type=int) 
@click.option("--out")
def infer(**kwargs):
    from .infer import threads_infer
    threads_infer(**kwargs)
    goodbye()

@main.command()
@click.option("--scaffold", required=True, help="Path to vcf containing phased scaffold of common variants")
@click.option("--argn", help="Path to reference ARG in .argn format")
@click.option("--ts", help="Path to reference ARG in .ts format")
@click.option("--unphased", required=True, help="Path to vcf containing the full target dataset (including scaffold variants)")
@click.option("--out", required=True, help="Path to phased output vcf")
def phase(**kwargs):
    from .phase import threads_phase
    threads_phase(**kwargs)
    goodbye()

@main.command()
@click.option("--threads", required=True, help="Path to an input .threads file")
@click.option("--argn", default=None, help="Path to an output .argn file")
@click.option("--tsz", default=None, help="Path to an output .tsz file")
@click.option("--add_mutations", is_flag=True, default=False, help="If passed, mutations are parsimoniously added to the output ARG. This may result in a high number of mutations if the --fit_to_data flag was not used.")
def convert(**kwargs):
    from .convert import threads_convert
    threads_convert(**kwargs)
    goodbye()

@main.command()
@click.option("--threads", required=True, help="Path to an input .threads file.")
@click.option("--out", required=True, help="Path to output.")
@click.option("--num_threads", type=int, help="Size of processor pool to process batches", default=None)
def allele_ages(**kwargs):
    from .allele_ages import estimate_allele_ages
    estimate_allele_ages(**kwargs)
    goodbye()

@main.command()
@click.option("--argn", help="Path to input .argn file")
@click.option("--out", help="Path to output .mut file")
@click.option("--maf", type=float, default=0.02, help="Do not store entries with MAF above this")
@click.option("--input", type=str, help="Path to bcf/vcf with genotypes to map with AC/AN fields")
@click.option("--region", type=str, help="Region in chr:start-end format (start and end inclusive)")
@click.option("--num_threads", type=int, help="Number of computational threads to request", default=1)
def map(**kwargs):
    from .map_mutations_to_arg import threads_map_mutations_to_arg
    threads_map_mutations_to_arg(**kwargs)
    goodbye()

@main.command()
@click.option("--panel", required=True, help="pgen array panel")
@click.option("--target", required=True, help="pgen array targets")
@click.option("--mut", required=True, help="pgen array targets")
@click.option("--map", required=True, help="Path to genotype map in SHAPEIT format")
@click.option("--mutation_rate", type=float, help="Per-site-per-generation SNP mutation rate", default=1.4e-8)
@click.option("--demography", required=True, help="Path to file containing demographic history")
@click.option("--out", help="Path to output .vcf file", default=None)
@click.option("--stdout", help="Redirect output to stdout (will disable logging)", is_flag=True)
@click.option("--region", required=True, type=str, help="Region in chr:start-end format (start and end inclusive)")
def impute(panel, target, map, mut, demography, out, stdout, region, mutation_rate=1.4e-8):
    # --stdout flag is mutually exclusive to --out flag. It is used only here to
    # confirm the user wants to redirect (potentially a lot of data) to stdout.
    # The Impute class does not use this variable, instead 'out' is just None.
    import sys
    if (stdout and out) or not (stdout or out):
        print("Either --out or --stdout must be specified", file=sys.stderr)
        exit(1)

    from .impute import Impute
    Impute(panel, target, map, mut, demography, out, region, mutation_rate)

    # Do not print anything in stdout mode, to keep output clean.
    if not stdout:
        goodbye()

@main.command()
@click.argument("threads", required=True)
def vcf(threads):
    """Convert THREADS to VCF format and print to stdout."""
    from .threads_to_vcf import threads_to_vcf
    threads_to_vcf(threads)

if __name__ == "__main__":
    main()
