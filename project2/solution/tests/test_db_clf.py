# HK, 2020-12-05


from unittest import TestCase

from tcga.utils import download

from idiva.clf.df import join, c5_df
from idiva.db.db_clf import classify
from idiva.download import download
from idiva.io import cache_df
from idiva.io.gz import open_maybe_gz
from idiva.io.vcf import ReadVCF

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


def case_control_maker():
    dfs = {}
    for group, value in URLS.items():
        with download(value).now.open(mode='rb') as fd:
            with open_maybe_gz(fd, mode='r') as fd:
                dfs[group] = c5_df(ReadVCF(fd))

    return join(case=dfs['case'], ctrl=dfs['ctrl'], on=['CHROM', 'POS', 'REF', 'ALT', 'ID'])


class TestDbClf(TestCase):
    def test_db_clf(self):
        with download(URLS['ctrl']).now.open(mode='rb') as ctrl, download(URLS['case']).now.open(mode='rb') as case:
            with open_maybe_gz(ctrl, mode='r') as ctrl, open_maybe_gz(case, mode='r') as case:
                result = classify(case=ReadVCF(case), ctrl=ReadVCF(ctrl),
                                  case_control=cache_df(name=("case_control"), key='c5_df',
                                                        df_maker=case_control_maker))
        self.assertTrue(len(result.df))
