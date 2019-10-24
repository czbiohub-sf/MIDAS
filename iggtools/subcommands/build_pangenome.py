import multiprocessing
from collections import defaultdict
import Bio.SeqIO
from iggtools.common.argparser import add_subcommand
from iggtools.common.utils import tsprint, InputStream, parse_table, decoded, retry, suppress_exceptions, command
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


@suppress_exceptions
def clean_genes_noexc(genome_id):
    return clean_genes(genome_id)


def clean_all_genes(species_genomes_ids, clean_func, num_procs):
    p = multiprocessing.Pool(num_procs)
    temp_files = p.map(clean_func, species_genomes_ids)
    failed = []
    cleaned = {}
    for genome_id, files in zip(species_genomes_ids, temp_files):
        if files:
            cleaned[genome_id] = files
        else:
            failed.append(genome_id)
    return cleaned, failed


def clean_all_genes_with_retry(species_genomes_ids, species_id):
    num_procs = (multiprocessing.cpu_count() + 1) // 2
    cleaned, failed = clean_all_genes(species_genomes_ids, clean_genes_noexc, num_procs)
    if failed:
        tsprint(f"Failed for {len(failed)} genomes from species {species_id}, including {', '.join(failed[:10])} even with individual genome retries, when running with {num_procs} subprocs.  Issuing a final wave of retries with 4x fewer subprocs.")
        # This time, a failure will raise an exception.  So there can be no failures.
        more_cleaned, failed_again = clean_all_genes(failed, clean_genes, (num_procs + 3) // 4)
        assert not failed_again
        tsprint(f"The wave of retries for {len(failed)} failed genomes from species {species_id} was a success.")
        cleaned.update(more_cleaned)
    tsprint(f"Successfully cleaned all {len(species_genomes_ids)} genomes from species {species_id}")
    return reordered_dict(cleaned, species_genomes_ids)


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

    cleaned = clean_all_genes_with_retry(species_genomes_ids, species_id)

    command("rm -f genes.ffn")

    files_in_chunk = 10  # keep "cat" commands short
    cleaned_files = list(cleaned.values())
    for chunk_start in range(0, len(cleaned_files), files_in_chunk):
        chunk_files = cleaned_files[chunk_start : chunk_start + files_in_chunk]
        command("cat " + " ".join(str(genes_file) for genes_file, _ in chunk_files) + " >> genes.ffn")
        command("cat " + " ".join(str(genes_len) for _, genes_len in chunk_files) + " >> genes.len")





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
