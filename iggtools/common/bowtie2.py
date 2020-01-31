#!/usr/bin/env python3
import os
from iggtools.common.utils import tsprint, num_physical_cores, command, split


def build_bowtie2_db(bt2_db_dir, bt2_db_name, downloaded_files, threads=1):
    """
    Build Bowtie2 database of representative genomes or centroid genes
    for the species present in the sample, e.g. repgenomes OR pangenomes
    """
    bt2_db_suffixes = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
    if all(os.path.exists(f"{bt2_db_dir}/{bt2_db_name}.{ext}") for ext in bt2_db_suffixes):
        tsprint("Skipping bowtie2-build as database files appear to exist.")
        return
    command(f"rm -f {bt2_db_dir}/{bt2_db_name}.fa")
    command(f"touch {bt2_db_dir}/{bt2_db_name}.fa")

    for files in split(downloaded_files.values(), 20):  # keep "cat" commands short
        command("cat " + " ".join(files) + f" >> {bt2_db_dir}/{bt2_db_name}.fa")

    command(f"bowtie2-build --threads {threads} {bt2_db_dir}/{bt2_db_name}.fa {bt2_db_dir}/{bt2_db_name} > {bt2_db_dir}/bowtie2-build.log")


def bowtie2_align(args, bt2_db_dir, bt2_db_name, sort_aln=False, threads=1):
    """
    Use Bowtie2 to map reads to specified representative genomes or
    collections of centroids genes for the pangenome flow.
    """

    if args.debug and os.path.exists(f"{bt2_db_dir}/{bt2_db_name}.bam"):
        tsprint(f"Skipping Bowtie2 alignment in debug mode as temporary data exists: {bt2_db_dir}/{bt2_db_name}.bam")
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
        bt2_command = f"bowtie2 --no-unal -x {bt2_db_dir}/{bt2_db_name} {max_reads} --{aln_mode} --{aln_speed} --threads {threads} -q {r1} {r2}"
        if sort_aln:
            command(f"set -o pipefail; {bt2_command} | \
                    samtools view --threads {threads} -b - | \
                    samtools sort --threads {threads} -o {bt2_db_dir}/{bt2_db_name}.bam")
        else:
            command(f"set -o pipefail; {bt2_command} | \
                    samtools view --threads {threads} -b - > {bt2_db_dir}/{bt2_db_name}.bam")
    except:
        tsprint(f"Bowtie2 align to {bt2_db_dir}/{bt2_db_name}.bam run into error")
        command(f"rm -f {bt2_db_dir}/{bt2_db_name}.bam")
        raise


def samtools_index(args, bt2_db_dir, bt2_db_name, threads):
    if args.debug and os.path.exists(f"{bt2_db_dir}/{bt2_db_name}.bam.bai"):
        tsprint(f"Skipping samtools index in debug mode as temporary data exists: {bt2_db_dir}/{bt2_db_name}.bam")
        return

    try:
        command(f"samtools index -@ {threads} {bt2_db_dir}/{bt2_db_name}.bam")
    except:
        command(f"rm -f {bt2_db_dir}/{bt2_db_name}.bam.bai")
        raise
