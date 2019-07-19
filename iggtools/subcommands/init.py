from iggtools.common import argparser
from iggtools.common.utils import tsprint, find_files, InputStream, OutputStream, parse_table
from iggtools.params import inputs, outputs


def init(args):
    """
    Input spec: https://github.com/czbiohub/iggtools/wiki#inputs
    Output spec: https://github.com/czbiohub/iggtools/wiki#target-layout-in-s3
    """

    msg = f"Building {outputs.genomes}."
    if find_files(outputs.genomes):
        if not args.force:
            tsprint(f"Destination {outputs.genomes} already exists.  Specify --force to overwrite.")
            return
        msg = f"Rebuilding {outputs.genomes}."
    tsprint(msg)

    id_remap = {}
    with InputStream(inputs.alt_species_ids) as ids:
        for row in parse_table(ids, ["alt_species_id", "species_id"]):
            new_id, old_id = row
            id_remap[old_id] = new_id


    seen_genomes, seen_species = set(), set()
    with OutputStream(outputs.genomes) as out:

        target_columns = ["genome", "species", "representative", "genome_is_representative"]
        out.write(("\t".join(target_columns) + "\n").encode('utf-8'))

        with InputStream(inputs.genomes2species) as g2s:
            for row in parse_table(g2s, ["MAG_code", "Species_id"]):
                genome, representative = row
                species = id_remap[representative]
                genome_is_representative = str(int(genome == representative))
                target_row = [genome, species, representative, genome_is_representative]
                out.write(("\t".join(target_row) + "\n").encode('utf-8'))
                seen_genomes.add(genome)
                seen_species.add(species)

    tsprint(f"Emitted {len(seen_genomes)} genomes and {len(seen_species)} species to {outputs.genomes}.")


def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand}.")
    assert args.subcommand == "init"
    init(args)


def register_argparser(singleton=[None]):  # pylint: disable=dangerous-default-value
    subparser = singleton[0]
    if not subparser:
        subparser = argparser.get().subparsers.add_parser('init', help=f"initialize target {outputs.genomes}")
        subparser.set_defaults(subcommand_main=main)
        argparser.add_shared_subcommand_args(subparser)
        singleton[0] = subparser


register_argparser()
