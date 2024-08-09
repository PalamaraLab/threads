# This file is part of the Threads software suite.
# Copyright (C) 2024 Threads Developers.
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
       `=`,'=_______________________________=`,'=/
         y== ||Thank you for using Threads|| y==/
       ,=,-<=‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾=,-<=`.
    ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
    """
    )


@click.group()
def main():
    pass


@main.command()
@click.option("--pgen", required=True, help="Path to input genotypes in pgen format.")
@click.option("--map_gz", default=None, help="Path to input genotype map with columns chromosome, snp, cM-position, bp-position.")
@click.option("--recombination_rate", default=1.3e-8, type=float, help="Genome-wide recombination rate. Ignored if a map is passed.")
@click.option("--demography", required=True, help="Path to input genotype.")
@click.option("--mode", required=True, type=click.Choice(['array', 'wgs']), default="wgs", help="Inference mode (wgs or array).")
@click.option("--query_interval", type=float, default=0.01, help="Hyperparameter for the preliminary haplotype matching in cM.")
@click.option("--match_group_interval", type=float, default=0.5, help="Hyperparameter for the preliminary haplotype matching in cM.")
@click.option("--mutation_rate", required=True, type=float, default=1.4e-8, help="Genome-wide mutation rate.")
@click.option("--num_threads", type=int, default=1, help="Number of computational threads to request.")
@click.option("--region", help="Region of genome for which ARG is output. The full genotype is still used for inference.")
@click.option("--max_sample_batch_size", help="Max number of LS processes run simultaneously per thread.", default=None, type=int) 
@click.option("--out")
def infer(pgen, map_gz, recombination_rate, demography, mutation_rate, query_interval, match_group_interval, mode, num_threads, region, max_sample_batch_size, out):
    from .infer import threads_infer
    threads_infer(pgen, map_gz, recombination_rate, demography, mutation_rate, query_interval, match_group_interval, mode, num_threads, region, max_sample_batch_size, out)
    goodbye()


@click.command()
@click.option("--scaffold", required=True, help="Path to vcf containing phased scaffold of common variants")
@click.option("--argn", help="Path to reference ARG in .argn format")
@click.option("--ts", help="Path to reference ARG in .ts format")
@click.option("--unphased", required=True, help="Path to vcf containing the full target dataset (including scaffold variants)")
@click.option("--out", required=True, help="Path to phased output vcf")
def phase(scaffold, argn, ts, unphased, out):
    from .phase import threads_phase
    threads_phase(scaffold, argn, ts, unphased, out)
    goodbye()


@main.command()
@click.option("--threads", required=True, help="Path to an input .threads file.")
@click.option("--argn", default=None, help="Path to an output .argn file.")
@click.option("--tsz", default=None, help="Path to an output .tsz file.")
@click.option("--max_n", default=None, help="How many samples to thread.", type=int)
@click.option("--random-seed", default=None, help="Seed for noise generation.", type=int)
@click.option("--verify", is_flag=True, show_default=True, default=False, help="Whether to use tskit to verify the ARG.")
def convert(threads, argn, tsz, max_n, random_seed, verify):
    from .convert import threads_convert
    threads_convert(threads, argn, tsz, max_n, random_seed, verify)
    goodbye()


@main.command()
def map():
    print("threads map implementation coming soon")
    goodbye()


@main.command()
@click.option("--panel", required=True, help="pgen array panel")
@click.option("--target", required=True, help="pgen array targets")
@click.option("--mut", required=True, help="pgen array targets")
@click.option("--map", required=True, help="Path to genotype map")
@click.option("--mutation_rate", type=float, help="Per-site-per-generation SNP mutation rate.", default=1.4e-8)
@click.option("--demography", required=True, help="Path to file containing demographic history.")
@click.option("--out", required=True, help="Path to output .vcf file.")
@click.option("--region", required=True, type=str)
def impute(panel, target, map, mut, demography, out, region, mutation_rate=1.4e-8):
    from .impute import Impute
    Impute(panel, target, map, mut, demography, out, region, mutation_rate)
    goodbye()


@click.command()
@click.argument("mode", type=click.Choice(["infer", "convert", "map"]))
@click.option("--argn", help="Path to .argn file with results")
@click.option("--out", help="Where to save results")
@click.option("--maf", type=float, default=0.02, help="Don't map stuff with MAF above this")
@click.option("--input", type=str, help="Path to bcf/vcf with genotypes to map. Most have AC/AN fields")
@click.option("--region", type=str, help="Of format chr:start-end (both inclusive)")
@click.option("--threads", type=int, help="Of format chr:start-end (both inclusive)", default=1)
def map_mutations_to_arg(argn, out, maf, input, region, threads):
    from .map_mutations_to_arg import threads_map_mutations_to_arg
    threads_map_mutations_to_arg(argn, out, maf, input, region, threads)
    goodbye()


if __name__ == "__main__":
    main()
