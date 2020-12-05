# HK, 2020-12-05


from idiva import log

import pandas
import io

from unittest import TestCase
from pathlib import Path
from tcga.utils import download
from idiva.io import open_maybe_gz
from idiva.db.db_clf import classify
from idiva.io.vcf import ReadVCF
from idiva.io import open_maybe_gz

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control_v2.vcf",
    'case': BASE / "case_processed_v2.vcf",
}


class TestDbClf(TestCase):
    def test_db_clf(self):
        with open(PATHS['ctrl'], mode='r') as control, open(PATHS['case'], mode='r') as case:
        # with open_maybe_gz(PATHS['ctrl']) as control, open_maybe_gz(PATHS['case']) as case:
            classify(case=ReadVCF(case), ctrl=ReadVCF(control))
