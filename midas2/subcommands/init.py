from midas2.common.argparser import add_subcommand
from midas2.common.utils import tsprint, find_files, InputStream, OutputStream, select_from_tsv
from midas2.params import inputs, outputs

def init(args):
    """
    Input spec: https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#inputs
    Output spec: https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB#target-layout-in-s3
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
        for row in select_from_tsv(ids, selected_columns=["alt_species_id", "species_id"]):
            new_id, old_id = row
            id_remap[old_id] = new_id

    seen_genomes, seen_species = set(), set()
    with OutputStream(outputs.genomes) as out:

        target_columns = ["genome", "species", "representative", "genome_is_representative"]
        out.write("\t".join(target_columns) + "\n")

        with InputStream(inputs.genomes2species) as g2s:
            for row in select_from_tsv(g2s, selected_columns=["MAG_code", "Species_id"]):
                genome, representative = row
                species = id_remap[representative]
                genome_is_representative = str(int(genome == representative))
                target_row = [genome, species, representative, genome_is_representative]
                out.write("\t".join(target_row) + "\n")
                seen_genomes.add(genome)
                seen_species.add(species)

    tsprint(f"Emitted {len(seen_genomes)} genomes and {len(seen_species)} species to {outputs.genomes}.")


def register_args(main_func):
    add_subcommand('init', main_func, help=f"initialize target {outputs.genomes}", epilog=init.__doc__)
    return main_func


@register_args
def main(args):
    tsprint(f"Executing midas2 subcommand {args.subcommand} with args {vars(args)}.")
    init(args)
