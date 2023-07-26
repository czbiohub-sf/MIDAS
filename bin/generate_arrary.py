import re
from midas2.models.midasdb import MIDAS_DB
from midas2.common.utils import OutputStream, command
midas_db = MIDAS_DB("/pollard/data/midas2-db/midas2db-wis-2023", "newdb", 16)
genomes = midas_db.uhgg.genomes
command("mkdir -p sge_jobs/step1_annotate")
n = 500
for i in range(0, 500):
    selected_genomes = set()
    for gid in genomes:
        gid_int = int(re.search(r'\d+$', gid).group())
        if gid_int % n == i:
            selected_genomes.add(gid)
    j = i + 1
    with OutputStream(f"sge_jobs/step1_annotate/batch_{j}") as stream:
        stream.write("\n".join(selected_genomes) + "\n")
