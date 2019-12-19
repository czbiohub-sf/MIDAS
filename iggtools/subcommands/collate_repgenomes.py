import os
from collections import defaultdict
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, parse_table, retry, command, multithreading_map, find_files, sorted_dict, upload, num_physical_cores, split
from iggtools.params import inputs, outputs


CONCURRENT_MARKER_GENES_IDENTIFY = num_physical_cores


def imported_genome_file(genome_id, species_id, component):
    return f"{outputs.cleaned_imports}/{species_id}/{genome_id}/{component}"


def input_marker_genes_file(genome_id, species_id, filename):
    # s3://{igg}/marker_genes/phyeco/temp/{SPECIES_ID}/{GENOME_ID}/{GENOME_ID}.{hmmsearch, markers.fa, markers.map}
    return f"{outputs.marker_genes}/temp/{species_id}/{genome_id}/{filename}"


def output_marker_genes():
    return f"{outputs.marker_genes}/{inputs.marker_set}.fa"


def read_toc(genomes_tsv, deep_sort=False):
    # Read in the table of contents.
    # We will centralize and improve this soon.
    species = defaultdict(dict)
    representatives = {}
    genomes = {}
    with InputStream(genomes_tsv) as table_of_contents:
        for row in parse_table(table_of_contents, ["genome", "species", "representative", "genome_is_representative"]):
            genome_id, species_id, representative_id, _ = row
            species[species_id][genome_id] = row
            representatives[species_id] = representative_id
            genomes[genome_id] = species_id
    if deep_sort:
        for sid in species.keys():
            species[sid] = sorted_dict(species[sid])
        species = sorted_dict(species)
    return species, representatives, genomes


@retry
def find_files_with_retry(f):
    return find_files(f)


def drop_lz4(filename):
    assert filename.endswith(".lz4")
    return filename[:-4]


# 1. Occasional failures in aws s3 cp require a retry.
@retry
def download_reference(ref_path_local_base):
    ref_path, local_base = ref_path_local_base
    local_path = os.path.basename(ref_path)
    local_path = drop_lz4(local_path)
    local_path = f"{local_base}/{local_path}"
    #command(f"rm -f {local_path}")
    #command(f"aws s3 cp --only-show-errors {ref_path} - | lz4 -dc > {local_path}")
    return local_path


def collate_repgenomes(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = os.path.basename(outputs.genomes)
    command(f"rm -f {local_toc}")
    command(f"aws s3 cp --only-show-errors {outputs.genomes} {local_toc}")

    species, representatives, _ = read_toc(local_toc)

    dest_file = output_marker_genes()
    msg = f"Collating marker genes sequences."
    if find_files_with_retry(dest_file):
        if not args.force:
            tsprint(f"Destination {dest_file} already exists.  Specify --force to overwrite.")
            return
        msg = msg.replace("Collating", "Recollating")
    tsprint(msg)

    slave_log = "collate_repgenomes.log"
    slave_subdir = f"collate_repgenomes"
    if not args.debug:
        command(f"rm -rf {slave_subdir}")
    if not os.path.isdir(slave_subdir):
        command(f"mkdir {slave_subdir}")
    slave_cmd = f"cd {slave_subdir}"
    with open(f"{slave_subdir}/{slave_log}", "w") as slog:
        slog.write(msg + "\n")
        slog.write(slave_cmd + "\n")
    command(slave_cmd)

    # Download
    ref_path_list = []
    for species_id in species.keys():
        representative_id = representatives[species_id]
        s3_marker_path = input_marker_genes_file(representative_id, species_id, f"{representative_id}.markers.fa.lz4")
        ref_path_list.append((s3_marker_path, slave_subdir))
    downloaded_markers = multithreading_map(download_reference, ref_path_list, num_threads=20)

    ## Collate
    command(slave_cmd)
    local_dest_file = os.path.basename(dest_file)
    downloaded_markers_basename = [os.path.basename(dm) for dm in downloaded_markers]
    for fna_files in split(downloaded_markers_basename, 20):  # keep "cat" commands short
        #command("cat " + " ".join(fna_files) + f" >> {local_dest_file}")
        command("cat " + " ".join(fna_files) + " >> phyeco.fa")

    ## Upload
    upload(f"{local_dest_file}", f"{dest_file}", check=False)

    ## Clean up
    if not args.debug:
        command(f"rm -rf {slave_subdir}", check=False)


def register_args(main_func):
    subparser = add_subcommand('collate_repgenomes', main_func, help='collate marker genes for repgresentative genomes')
    subparser.add_argument('--species',
                           dest='species',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning import genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    collate_repgenomes(args)
