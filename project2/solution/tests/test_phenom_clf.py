# HK, 2020-12-07


from unittest import TestCase

from idiva.clf.phenomenet_classifier import phenomenet_classifier
from idiva.clf.phenomenet_clf_basic import phenomenet_classifier_basic
from idiva.io.vcf import ReadVCF
from idiva import log

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestPhenomClf(TestCase):
    def test_phenom_clf(self):
        with ReadVCF.open(URLS['case']) as case, ReadVCF.open(URLS['ctrl']) as ctrl:
            result = phenomenet_classifier(case=case, ctrl=ctrl)
        self.assertTrue(len(result.df))
        log.info('passed!')

    def test_phenom_basic(self):
        with ReadVCF.open(URLS['case']) as case:
            result = phenomenet_classifier_basic(case=case, ctrl=None)
        self.assertTrue(len(result.df))
        log.info('passed!')
