# HK, 2020-12-05


from unittest import TestCase

from idiva.db.db_clf import db_classifier
from idiva.io.vcf import ReadVCF

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestDbClf(TestCase):
    def test_db_clf(self):
        with ReadVCF.open(URLS['case']) as case:
            result = db_classifier(case=case, ctrl=None)
        self.assertTrue(len(result.df))
