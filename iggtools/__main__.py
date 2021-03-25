#!/usr/bin/env python3
#
# Import every subcommand you wish to expose.  Omitting a subcommand here will hide it.
#
# Pylint thinks these modules are "unused" but importing them has the side effect
# of registering their respective subcommands into the argument parser.
#
# Each subcommand module has its own main() function.  They are totally independent,
# with their own set of command line arguments and subcommand help text -- aside from
# the shared arguments and shared help text defined in iggtools.common.argparser.
#
from iggtools.subcommands import aws_batch_init, aws_batch_submit, init, \
    import_uhgg, annotate_genome, build_pangenome, \
    infer_markers, build_markerdb, generate_dbmisc, \
    midas_run_species, midas_run_genes, midas_run_snps, \
    midas_merge_species, midas_merge_snps, midas_merge_genes, \
    build_bowtie2_indexes # pylint: disable=unused-import

from iggtools.common.argparser import parse_args


def main():
    subcommand_main, subcommand_args = parse_args()
    return subcommand_main(subcommand_args)


if __name__ == "__main__":
    main()
