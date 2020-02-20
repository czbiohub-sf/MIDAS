# Base class for Sample: map of directory structure of a sample
import os
from collections import defaultdict
from iggtools.common.utils import InputStream, select_from_tsv
from iggtools.params.schemas import fetch_schema_by_dbtype


class Sample: # pylint: disable=too-few-public-methods

    def __init__(self, sample_name, midas_outdir, dbtype=None):
        self.sample_name = sample_name
        self.midas_outdir = midas_outdir
        assert os.path.exists(midas_outdir), f"Provided MIDAS output {midas_outdir} for {sample_name} in sample_list is invalid"

        self.data_dir = os.path.join(self.midas_outdir, dbtype)
        assert os.path.exists(self.data_dir), f"Missing MIDAS {dbtype} directiory for {self.data_dir} for {sample_name}"

        self.profile = {}
        if dbtype is not None:
            self.profile = self.fetch_profile(dbtype)

    def fetch_profile(self, dbtype):
        summary_path = f"{self.data_dir}/summary.txt"
        assert os.path.exists(summary_path), f"Missing MIDAS {summary_path} for {self.sample_name}"

        schema = fetch_schema_by_dbtype(dbtype)
        profile = defaultdict()
        with InputStream(self.profile_path) as stream:
            for info in select_from_tsv(stream, selected_columns=schema, result_structure=dict):
                profile[info["species_id"]] = info
        return profile
