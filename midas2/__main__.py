#!/usr/bin/env python3
#
# Import every subcommand you wish to expose.  Omitting a subcommand here will hide it.
#
# Pylint thinks these modules are "unused" but importing them has the side effect
# of registering their respective subcommands into the argument parser.
#
# Each subcommand module has its own main() function.  They are totally independent,
# with their own set of command line arguments and subcommand help text -- aside from
# the shared arguments and shared help text defined in midas2.common.argparser.
#
from midas2.subcommands import aws_batch_init, aws_batch_submit, init, \
    import_uhgg, annotate_genome, build_pangenome, \
    infer_markers, build_midasdb, database, \
    run_species, run_genes, run_snps, \
    merge_species, merge_snps, merge_genes, \
    build_bowtie2db, compute_chunks # pylint: disable=unused-import

from midas2.common.argparser import parse_args


def main():
    subcommand_main, subcommand_args = parse_args()
    return subcommand_main(subcommand_args)


if __name__ == "__main__":
    main()
