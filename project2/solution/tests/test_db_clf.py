# HK, 2020-12-05


from unittest import TestCase

from tcga.utils import download

from idiva.db.db_clf import db_classifier
from idiva.download import download
from idiva.io.gz import open_maybe_gz
from idiva.io.vcf import ReadVCF
from idiva.io import cache_df
from idiva.clf.df import c5_df
import pandas as pd

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestDbClf(TestCase):
    def test_db_clf(self):
        with ReadVCF.open(URLS['case']) as case:
            result = db_classifier(case=case, ctrl=None)
        self.assertTrue(len(result.df))
