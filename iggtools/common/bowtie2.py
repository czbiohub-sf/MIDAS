#!/usr/bin/env python3
import os
from iggtools.common.utils import tsprint, num_physical_cores, command, split



def build_bowtie2_db(bt2_db_dir, bt2_db_name, downloaded_files):
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

    command(f"bowtie2-build --threads {num_physical_cores} {bt2_db_dir}/{bt2_db_name}.fa {bt2_db_dir}/bt2_db_name > {bt2_db_dir}/bowtie2-build.log")
