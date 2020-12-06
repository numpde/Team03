# HK, 2020-12-05


from unittest import TestCase

from tcga.utils import download

from idiva.db.db_clf import db_classifier
from idiva.download import download
from idiva.io.gz import open_maybe_gz
from idiva.io.vcf import ReadVCF

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestDbClf(TestCase):
    def test_db_clf(self):
        with download(URLS['ctrl']).now.open(mode='rb') as ctrl, download(URLS['case']).now.open(mode='rb') as case:
            with open_maybe_gz(ctrl, mode='r') as ctrl, open_maybe_gz(case, mode='r') as case:
                result = db_classifier(case=ReadVCF(case), ctrl=ReadVCF(ctrl))
        self.assertTrue(len(result.df))
