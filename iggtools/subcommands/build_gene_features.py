import os
import sys
from multiprocessing import Semaphore
from iggtools.common.argparser import add_subcommand, SUPPRESS
from iggtools.common.utils import tsprint, retry, command, multithreading_map, find_files, upload, pythonpath, upload_star, download_reference, InputStream, OutputStream
from iggtools.models.uhgg import UHGG
from iggtools.params import outputs
from iggtools.subcommands.import_uhgg import decode_genomes_arg
import gffutils
import Bio.SeqIO


# Up to this many concurrent species builds.
CONCURRENT_GENOME_BUILDS = Semaphore(20)


def annotations_file(genome_id, species_id, filename):
    # s3://microbiome-igg/2.0/prodigal/GUT_GENOMEDDDDDD.{fna, faa, gff, log}
    return f"{outputs.annotations}/{species_id}/{genome_id}/{filename}"


@retry
def find_files_with_retry(f):
    return find_files(f)


## Where should I put this read in sequence file depends on the usage, let's circle back this on Thursday
#fasta_file = "/Users/chunyu.zhao/Desktop/20210104_test/GUT_GENOME000001/GUT_GENOME000001.ffn"
def read_gene_sequence(fasta_file, species_id):
    """ Scan the genome file to get contig_id and contig_seq as ref_seq """
    contigs = {}
    with InputStream(fasta_file) as file:
        for rec in Bio.SeqIO.parse(file, 'fasta'):
            contigs[rec.id] = {
                "species_id": species_id,
                "contig_id": rec.id,
                "contig_len": len(rec.seq),
                "contig_seq": str(rec.seq),
            }
    return contigs


#gff3_file = "/Users/chunyu.zhao/Desktop/20210104_test/GUT_GENOME000001/GUT_GENOME000001.gff"
def reformat_gene_features(gff3_file, genes_file):
    # Convert GFF3 features into format desired by MIDAS
    db = gffutils.create_db(gff3_file, f"{gff3_file}.db")

    with OutputStream(genes_file) as stream:
        stream.write("\t".join(["gene_id", "contig_id", "start", "end", "strand", "gene_type"]) + "\n")
        for feature in db.all_features():
            if feature.source == "prokka":
                continue
            seqid = feature.seqid.replace("gnl|Prokka|", "")
            start = feature.start - 1
            stop = feature.stop
            strand = feature.strand
            gene_id = feature.attributes['ID'][0]
            locus_tag = feature.attributes['locus_tag'][0]
            assert gene_id == locus_tag
            gene_type = feature.featuretype
            stream.write("\t".join([gene_id, seqid, str(start), str(stop), strand, gene_type]) + "\n")
    return True


def build_gene_features(args):
    if args.zzz_slave_toc:
        build_gene_features_slave(args)
    else:
        build_gene_features_master(args)


def build_gene_features_master(args):

    # Fetch table of contents from s3.
    # This will be read separately by each species build subcommand, so we make a local copy.
    local_toc = download_reference(outputs.genomes)
    db = UHGG(local_toc)
    species_for_genome = db.genomes


    def genome_work(genome_id):
        assert genome_id in species_for_genome, f"Genome {genome_id} is not in the database."
        species_id = species_for_genome[genome_id]

        # The species build will upload this file last, after everything else is successfully uploaded.
        # Therefore, if this file exists in s3, there is no need to redo the species build.
        dest_file = annotations_file(genome_id, species_id, f"{genome_id}.genes.lz4")
        msg = f"Builing gene features for genome {genome_id} from species {species_id}."
        if find_files_with_retry(dest_file):
            if not args.force:
                tsprint(f"Destination {dest_file} for genome {genome_id} gene features already exists.  Specify --force to overwrite.")
                return
            msg = msg.replace("Importing", "Reimporting")


        with CONCURRENT_GENOME_BUILDS:
            tsprint(msg)
            slave_log = "build_gene_features.log"
            slave_subdir = f"{species_id}__{genome_id}"
            if not args.debug:
                command(f"rm -rf {slave_subdir}")
            if not os.path.isdir(slave_subdir):
                command(f"mkdir {slave_subdir}")

            # Recurisve call via subcommand.  Use subdir, redirect logs.
            slave_cmd = f"cd {slave_subdir}; PYTHONPATH={pythonpath()} {sys.executable} -m iggtools build_gene_features --genome {genome_id} --zzz_slave_mode --zzz_slave_toc {os.path.abspath(local_toc)} {'--debug' if args.debug else ''} &>> {slave_log}"
            with open(f"{slave_subdir}/{slave_log}", "w") as slog:
                slog.write(msg + "\n")
                slog.write(slave_cmd + "\n")
            try:
                command(slave_cmd)
            finally:
                # Cleanup should not raise exceptions of its own, so as not to interfere with any
                # prior exceptions that may be more informative.  Hence check=False.
                upload(f"{slave_subdir}/{slave_log}", annotations_file(genome_id, species_id, slave_log + ".lz4"), check=False)
                if not args.debug:
                    command(f"rm -rf {slave_subdir}", check=False)

    genome_id_list = decode_genomes_arg(args, species_for_genome)
    multithreading_map(genome_work, genome_id_list, num_threads=20)


def build_gene_features_slave(args):
    """
    https://github.com/czbiohub/iggtools/wiki
    """

    violation = "Please do not call build_gene_features_slave directly.  Violation"
    assert args.zzz_slave_mode, f"{violation}:  Missing --zzz_slave_mode arg."
    assert os.path.isfile(args.zzz_slave_toc), f"{violation}: File does not exist: {args.zzz_slave_toc}"

    db = UHGG(args.zzz_slave_toc)
    species_for_genome = db.genomes

    genome_id = args.genomes
    species_id = species_for_genome[genome_id]

    annot_file = f"{genome_id}.gff"
    download_reference(annotations_file(genome_id, species_id, f"{annot_file}.lz4"))

    out_file = f"{genome_id}.genes"
    dest_file = annotations_file(genome_id, species_id, f"{out_file}.lz4")
    #command(f"aws s3 rm --recursive {dest_file.rsplit('/', 1)[0]}")

    assert reformat_gene_features(annot_file, out_file)
    upload(out_file, dest_file)


def register_args(main_func):
    subparser = add_subcommand('build_gene_features', main_func, help='convert GFF into desired format')
    subparser.add_argument('--genomes',
                           dest='genomes',
                           required=False,
                           help="genome[,genome...] to import;  alternatively, slice in format idx:modulus, e.g. 1:30, meaning annotate genomes whose ids are 1 mod 30; or, the special keyword 'all' meaning all genomes")
    subparser.add_argument('--zzz_slave_toc',
                           dest='zzz_slave_toc',
                           required=False,
                           help=SUPPRESS) # "reserved to pass table of contents from master to slave"
    return main_func


@register_args
def main(args):
    tsprint(f"Executing iggtools subcommand {args.subcommand} with args {vars(args)}.")
    build_gene_features(args)
