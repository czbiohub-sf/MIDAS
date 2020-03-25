#!/usr/bin/env python3
import os
from iggtools.common.utils import tsprint, num_physical_cores, command, split


def bowtie2_index_exists(bt2_db_dir, bt2_db_name):
    bt2_db_suffixes = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
    if all(os.path.exists(f"{bt2_db_dir}/{bt2_db_name}.{ext}") for ext in bt2_db_suffixes):
        tsprint(f"Bowtie2 index {bt2_db_dir}/{bt2_db_name} exist.")
        return True
    return False


def build_bowtie2_db(bt2_db_dir, bt2_db_name, downloaded_files, cat=True):
    """ Build Bowtie2 database for the collections of fasta files """

    bt2_db_prefix = f"{bt2_db_dir}/{bt2_db_name}"
    if not bowtie2_index_exists(bt2_db_dir, bt2_db_name):
        if cat:
            command(f"rm -f {bt2_db_dir}/{bt2_db_name}.fa")
            command(f"touch {bt2_db_dir}/{bt2_db_name}.fa")
            for files in split(downloaded_files, 20):  # keep "cat" commands short
                command("cat " + " ".join(files) + f" >> {bt2_db_dir}/{bt2_db_name}.fa")
        else:
            command(f"ln -rs {downloaded_files} {bt2_db_dir}/{bt2_db_name}.fa")

        try:
            command(f"bowtie2-build --threads {num_physical_cores} {bt2_db_prefix}.fa {bt2_db_prefix} > {bt2_db_dir}/bowtie2-build.log")
        except:
            tsprint(f"Bowtie2 index {bt2_db_prefix} run into error")
            command(f"rm -f {bamfile_prefix}.1.bt2")
            raise

    return bt2_db_prefix


def bowtie2_align(bt2_db_dir, bt2_db_name, bamfile_path, args):
    """ Use Bowtie2 to map reads to prebuilt bowtie2 database """

    bt2_db_prefix = f"{bt2_db_dir}/{bt2_db_name}"
    bamfile_prefix = f"{bt2_db_dir}/{bt2_db_name}"

    if os.path.exists(bamfile_path):
        tsprint(f"Use existing per-sample bowtie2 db {bamfile_path}")
        return

    # Construct bowtie2 align input arguments
    max_reads = f"-u {args.max_reads}" if args.max_reads else ""
    aln_mode = "local" if args.aln_mode == "local" else "end-to-end"
    aln_speed = args.aln_speed if aln_mode == "end_to_end" else args.aln_speed + "-local"
    r2 = ""
    if args.r2:
        r1 = f"-1 {args.r1}"
        r2 = f"-2 {args.r2}"
    elif args.aln_interleaved:
        r1 = f"--interleaved {args.r1}"
    else:
        r1 = f"-U {args.r1}"

    try:
        bt2_command = f"bowtie2 --no-unal -x {bt2_db_prefix} {max_reads} --{aln_mode} --{aln_speed} --threads {num_physical_cores} -q {r1} {r2}"
        if args.aln_sort:
            command(f"set -o pipefail; {bt2_command} | \
                    samtools view --threads {num_physical_cores} -b - | \
                    samtools sort --threads {num_physical_cores} -o {bamfile_path}")
        else:
            command(f"set -o pipefail; {bt2_command} | \
                    samtools view --threads {num_physical_cores} -b - > {bamfile_path}")
    except:
        tsprint(f"Bowtie2 align to {bamfile_path} run into error")
        command(f"rm -f {bamfile_path}")
        raise


def samtools_index(bamfile_path, debug):

    if debug and os.path.exists(f"{bamfile_path}.bai"):
        tsprint(f"Skipping samtools index in debug mode as temporary data exists: {bamfile_path}.bai")
        return

    try:
        command(f"samtools index -@ {num_physical_cores} {bamfile_path}")
    except:
        tsprint(f"Samtools index {bamfile_path} run into error")
        command(f"rm -f {bamfile_path}.bai")
        raise
