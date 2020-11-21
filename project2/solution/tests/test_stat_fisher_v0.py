# RA, 2020-11-19

import pandas as pd
import numpy as np

from pathlib import Path
from unittest import TestCase
from warnings import filterwarnings

BASE = (Path(__file__).parent) / "data_for_tests"

PATHS_LARGE_HEAD = {
    'ctrl': BASE / "large_head/control_v2.vcf",
    'case': BASE / "large_head/case_processed_v2.vcf",
}

from tcga.utils import download

download_cache = (Path(__file__).parent.parent.parent / "input/download_cache")
assert download_cache.is_dir()
download = download.to(abs_path=download_cache)

URLS_LARGE = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}

MY_SPACE = BASE / F"my_space/{Path(__file__).stem}"
MY_SPACE.mkdir(exist_ok=True)

# RuntimeWarning: numpy.ufunc size changed
# https://github.com/numpy/numpy/issues/14920#issuecomment-554672523
filterwarnings("ignore", message="numpy.ufunc size changed")

from idiva.stat import v0_fisher
from scipy.stats import fisher_exact


class TestChi2(TestCase):
    def test_fisher_oneunit(self):
        (a1, a2) = ([1, 10, 30], [30, 10, 1])
        (c1, c2) = (["ALT0_case", "ALT1_case", "ALT2_case"], ["ALT0_ctrl", "ALT1_ctrl", "ALT2_ctrl"])
        df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])

        (oddrat0, pvalue0) = fisher_exact([[a1[0], a1[1] + a1[2]], [a2[0], a2[1] + a2[2]]])
        (oddrat1, pvalue1) = fisher_exact([[a1[1], a1[0] + a1[2]], [a2[1], a2[0] + a2[2]]])
        (oddrat2, pvalue2) = fisher_exact([[a1[2], a1[0] + a1[1]], [a2[2], a2[0] + a2[1]]])

        assert oddrat0 < 1e-2
        assert oddrat1 == 1
        assert oddrat2 > 1e+2

        pvalue_can = v0_fisher(df).to_numpy().squeeze().tolist()
        oddrat_can = v0_fisher(df, oddsratio=True).to_numpy().squeeze().tolist()

        pvalue_ref = [pvalue0, pvalue1, pvalue2]
        oddrat_ref = [oddrat0, oddrat1, oddrat2]

        r = (lambda L: list(np.round(x, 3) for x in L))

        self.assertListEqual(r(pvalue_can), r(pvalue_ref))
        self.assertListEqual(r(oddrat_can), r(oddrat_ref))


    # def test_chi2_na(self):
    #     (a1, a2) = ([5, 11, 20], [10, 11, pd.NA])
    #     (c1, c2) = (["A1", "B1", "C1"], ["A2", "B2", "C2"])
    #     df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])
    #     cdt = chi2_test(df, (c1, c2)).squeeze()
    #     self.assertTrue(pd.isna(cdt))

    def test_fisher_large_head(self):
        from idiva.io import ReadVCF, open_maybe_gz
        from idiva.clf.df import v0_df, join

        dfs = {}

        for (k, file) in PATHS_LARGE_HEAD.items():
            with open_maybe_gz(file, mode='r') as fd:
                dfs[k] = v0_df(ReadVCF(fd))

        df = join(case=dfs['case'], ctrl=dfs['ctrl'])

        pvalues = v0_fisher(df)
        print(pvalues.to_markdown())


    # def test_chi2_large(self):
    #     from idiva.io import ReadVCF, open_maybe_gz
    #     from idiva.clf.df import v0_df, join
    #
    #     dfs = {}
    #
    #     for k in URLS_LARGE:
    #         with download(URLS_LARGE[k]).now.open(mode='rb') as fd:
    #             with open_maybe_gz(fd, mode='r') as fd:
    #                 dfs[k] = v0_df(ReadVCF(fd))
    #
    #     df = join(case=dfs['case'], ctrl=dfs['ctrl'])
    #
    #     cols = tuple([F"ALT{n}_{kind}" for n in range(3)] for kind in ['case', 'ctrl'])
    #
    #     p = chi2_test(df[cols[0] + cols[1]], cols, add=1)
