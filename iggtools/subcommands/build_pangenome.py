import multiprocessing
from collections import defaultdict
from operator import itemgetter
import Bio.SeqIO
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, parse_table, decoded, retry, command, portion, multi_map
from iggtools.params import outputs


def pangenome_file(representative_id, component):
    # s3://microbiome-igg/2.0/pangenomes/GUT_GENOMEDDDDDD/{genes.ffn, centroids.ffn, gene_info.txt}
    return f"{outputs.pangenomes}/{representative_id}/{component}"


def annotations_file(genome_id, component_extension):
    # s3://microbiome-igg/2.0/prodigal/GUT_GENOMEDDDDDD.{fna, faa, gff, log}
    return f"{outputs.annotations}/{genome_id}.{component_extension}"


def sorted_dict(d):
    return {k: d[k] for k in sorted(d.keys())}


def reordered_dict(d, key_order):
    return {k: d[k] for k in key_order}


def read_toc(genomes_tsv, deep_sort=False):
    # Read in the table of contents.
    # We will centralize and improve this soon.
    species = defaultdict(dict)
    representatives = {}
    with InputStream(genomes_tsv) as table_of_contents:
        for row in parse_table(table_of_contents, ["genome", "species", "representative", "genome_is_representative"]):
            genome_id, species_id, representative_id, _ = row
            species[species_id][genome_id] = row
            representatives[species_id] = representative_id
    if deep_sort:
        for sid in species.keys():
            species[sid] = sorted_dict(species[sid])
        species = sorted_dict(species)
    return species, representatives


# aws s3e cp may occasionally fail, so should be retried once
# for really large numbers of genomes, we can also do a second pass of retries for the failures
@retry
def clean_genes(genome_id):
    input_annotations = annotations_file(genome_id, "fna.lz4")
    output_genes = f"{genome_id}.genes.fna"
    output_info = f"{genome_id}.genes.len"

    with open(output_genes, 'w') as o_genes, \
         open(output_info, 'w') as o_info, \
         InputStream(input_annotations, check_path=False) as genes:
        for rec in Bio.SeqIO.parse(decoded(genes), 'fasta'):
            gene_id = rec.id
            gene_seq = str(rec.seq).upper()
            gene_len = len(gene_seq)
            if gene_len == 0 or gene_id == '' or gene_id == '|':
                # Documentation for why we ignore these gene_ids should be added to
                # https://github.com/czbiohub/iggtools/wiki#pan-genomes
                # Also, we should probably count these and report stats.
                pass
            else:
                o_genes.write(f">{gene_id}\n{gene_seq}\n")
                o_info.write(f"{gene_id}\t{genome_id}\t{gene_len}\n")

    return output_genes, output_info


def build_pangenome(args):
    """
    Input spec:  https://github.com/czbiohub/iggtools/wiki#gene-annotations
    Output spec: https://github.com/czbiohub/iggtools/wiki#pan-genomes
    """

    species, representatives = read_toc(outputs.genomes)

    species_id = args.species
    assert species_id in species, f"Species {species_id} is not in the database."

    species_genomes = species[species_id]
    representative_id = representatives[species_id]
    species_genomes_ids = species_genomes.keys()
    tsprint(f"There are {len(species_genomes)} genomes for species {species_id} with representative genome {representative_id}.")

    num_procs = (multiprocessing.cpu_count() + 1) // 2
    cleaned = multi_map(clean_genes, species_genomes_ids, num_procs)

    command("rm -f genes.ffn")
    command("rm -f genes.len")

    for temp_files in portion(cleaned.values(), 20):  # keep "cat" commands short
        fna_files = map(itemgetter(0), temp_files)
        len_files = map(itemgetter(1), temp_files)
        command("cat " + " ".join(fna_files) + " >> genes.ffn")
        command("cat " + " ".join(len_files) + " >> genes.len")

    num_procs = (multiprocessing.cpu_count() + 1) // 2
    percents = [99, 90, 85, 80, 75]
    clusters = {pid: f"centroids.{pid}.ffn" for pid in percents}
    genes = "genes.ffn"
    for pid in percents:
        # Clusters "genes" into "clusters[pid]"
        # After the first iteration, the rest could run in parallel
        vsearch = f"vsearch -cluster_fast {genes} -id {pid/100.0} -threads {num_procs} -centroids {clusters[pid]} -uc uclust.{pid}.txt"
        genes = clusters[99] # fixed 99




def register_args(main_func):
    subparser = add_subcommand('build_pangenome', main_func, help='build pangenome for given species')
    subparser.add_argument('-s',
                           '--species',
                           dest='species',
                           required=True,
                           help="species whose pangenome to build")
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    build_pangenome(args)
