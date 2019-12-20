import os
from collections import defaultdict
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, parse_table, retry, command, multithreading_map, find_files, sorted_dict, upload, num_physical_cores, split, upload_star
from iggtools.params import inputs, outputs


CONCURRENT_MARKER_GENES_DOWNLOAD = num_physical_cores


def input_marker_genes_file(genome_id, species_id, filename):
    # s3://{igg}/marker_genes/phyeco/temp/{SPECIES_ID}/{GENOME_ID}/{GENOME_ID}.{hmmsearch, markers.fa, markers.map}
    return f"{outputs.marker_genes}/temp/{species_id}/{genome_id}/{filename}"


def output_all_rep_marker_genes(component):
    return f"{outputs.marker_genes}/{inputs.marker_set}.{component}"


def destpath(local_file):
    return f"{outputs.marker_genes}/{local_file}.lz4"


def check_dest_files(dest_file, msg, force):
    if find_files_with_retry(dest_file):
        if not force:
            tsprint(f"Destination {dest_file} already exists.  Specify --force to overwrite.")
            return msg
        msg = msg.replace(msg.split(" ")[0], "Re-" + msg.split(" ")[0])
    return msg


def read_toc(genomes_tsv, deep_sort=False):
    # Read in the table of contents.
    # We will centralize and improve this soon.
    species = defaultdict(dict)
    representatives = {}
    genomes = {}
    with InputStream(genomes_tsv) as table_of_contents:
        for row in parse_table(table_of_contents, ["genome", "species", "representative", "genome_is_representative"]):
            genome_id, species_id, rep_id, _ = row
            species[species_id][genome_id] = row
            representatives[species_id] = rep_id
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
    local_path = drop_lz4(os.path.basename(ref_path))
    local_path = f"{local_base}/{local_path}"
    command(f"rm -f {local_path}")
    command(f"aws s3 cp --only-show-errors {ref_path} - | lz4 -dc > {local_path}")
    return local_path


def collate_repgenomes(args):

    ## change the subcommand:collate_repgenomes_marker something like thatself.
    ## also need  to double check the phyeco.fa is the same with all the GUT genomes fa

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = os.path.basename(outputs.genomes)
    command(f"rm -f {local_toc}")
    command(f"aws s3 cp --only-show-errors {outputs.genomes} {local_toc}")

    species, representatives, _ = read_toc(local_toc)

    collate_log = "collate_repgenomes.log"
    collate_subdir = f"collate_repgenomes"

    last_output = destpath(collate_log)
    msg = f"Collating marker genes sequences."
    msg = check_dest_files(last_output, msg, args.force)

    tsprint(msg)
    with open(f"{collate_subdir}/{collate_log}", "w") as slog:
        slog.write(msg + "\n")

    if not args.debug:
        command(f"rm -rf {collate_subdir}")
    if not os.path.isdir(collate_subdir):
        command(f"mkdir {collate_subdir}")

    # Download
    download_seq_tasks = []
    download_map_tasks = []
    for species_id in species.keys():
        rep_id = representatives[species_id]
        remote_path_seq = input_marker_genes_file(rep_id, species_id, f"{rep_id}.markers.fa.lz4")
        remote_path_map = input_marker_genes_file(rep_id, species_id, f"{rep_id}.markers.map.lz4")
        download_seq_tasks.append((remote_path_seq, collate_subdir))
        download_map_tasks.append((remote_path_map, collate_subdir))
    downloaded_marker_seqs = multithreading_map(download_reference, download_seq_tasks, num_threads=CONCURRENT_MARKER_GENES_DOWNLOAD)
    downloaded_marker_maps = multithreading_map(download_reference, download_map_tasks, num_threads=CONCURRENT_MARKER_GENES_DOWNLOAD)

    ## Collate
    collated_rep_marker_seqs = output_all_rep_marker_genes("fa")
    collated_genes = os.path.basename(collated_rep_marker_seqs)
    for marker_fa_files in split(downloaded_marker_seqs, 20):
        command("cat " + " ".join(marker_fa_files) + f" >> {collate_subdir}/{collated_genes}")

    collated_rep_marker_maps = output_all_rep_marker_genes("map")
    collated_maps = os.path.basename(collated_rep_marker_maps)
    for marker_map_files in split(downloaded_marker_maps, 20):
        command("cat " + " ".join(marker_map_files) + f" >> {collate_subdir}/{collated_maps}")

    ## Index
    cmd_index = f"cd {collate_subdir}; hs-blastn index {collated_genes} &>> {collate_log}"
    with open(f"{collate_subdir}/{collate_log}", "a") as slog:
        slog.write(cmd_index + "\n")
    command(cmd_index)
    index_suffix = ["fa", "map", "fa.bwt", "fa.header", "fa.sa", "fa.sequence"]
    output_files = [f"{collate_subdir}/{inputs.marker_set}.{isuffix}" for isuffix in index_suffix]

    ## Upload
    upload_tasks = []
    for o in output_files:
        upload_tasks.append((o, destpath(os.path.basename(o))))
    multithreading_map(upload_star, upload_tasks)

    # Upload the log file in the last
    upload(f"{collate_subdir}/{collate_log}", destpath(collate_log), check=False)

    ## Clean up
    if not args.debug:
        command(f"rm -rf {collate_subdir}", check=False)


def register_args(main_func):
    subparser = add_subcommand('collate_repgenomes', main_func, help='collate marker genes for repgresentative genomes')
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    collate_repgenomes(args)
