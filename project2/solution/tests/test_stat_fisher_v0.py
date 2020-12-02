# RA, 2020-11-19

from idiva import log

import pandas as pd
import numpy as np
import io

from pathlib import Path
from unittest import TestCase

BASE = (Path(__file__).parent) / "data_for_tests"

PATHS_LARGE_HEAD = {
    'ctrl': BASE / "large_head/control_v2.vcf",
    'case': BASE / "large_head/case_processed_v2.vcf",
}

PATHS_LARGE_FULL = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}

MY_SPACE = BASE / F"my_space/{Path(__file__).stem}"
MY_SPACE.mkdir(exist_ok=True)

from idiva.stat import v0_fisher, fisher_scipy


class TestChi2(TestCase):
    def test_fisher_oneunit(self):
        (a1, a2) = ([1, 10, 20], [20, 10, 1])
        (c1, c2) = (["ALT0_case", "ALT1_case", "ALT2_case"], ["ALT0_ctrl", "ALT1_ctrl", "ALT2_ctrl"])
        df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])

        (oddrat0, pvalue0) = fisher_scipy([[a1[0], a1[1] + a1[2]], [a2[0], a2[1] + a2[2]]])
        (oddrat1, pvalue1) = fisher_scipy([[a1[1], a1[0] + a1[2]], [a2[1], a2[0] + a2[2]]])
        (oddrat2, pvalue2) = fisher_scipy([[a1[2], a1[0] + a1[1]], [a2[2], a2[0] + a2[1]]])

        assert oddrat0 < 1e-1
        assert oddrat1 == 1
        assert oddrat2 > 1e+1

        r = (lambda L: list(F"{x:.3e}" for x in L))

        pvalue_can = v0_fisher(df).to_numpy().squeeze().tolist()
        pvalue_ref = [pvalue0, pvalue1, pvalue2]

        self.assertListEqual(r(pvalue_can), r(pvalue_ref))

        with self.assertRaises(NotImplementedError):
            oddrat_can = v0_fisher(df, oddsratio=True).to_numpy().squeeze().tolist()
            oddrat_ref = [oddrat0, oddrat1, oddrat2]
            self.assertListEqual(r(oddrat_can), r(oddrat_ref))

    def test_fisher_large_head(self):
        from idiva.io import ReadVCF, open_maybe_gz

        from idiva.clf.df import v0_df
        from idiva.clf.df import join

        dfs = {}

        for (k, file) in PATHS_LARGE_HEAD.items():
            with open_maybe_gz(file, mode='r') as fd:
                assert isinstance(fd, io.TextIOBase)
                dfs[k] = v0_df(ReadVCF(fd))

        df = join(case=dfs['case'], ctrl=dfs['ctrl'])

        pvalues = v0_fisher(df)
        print(pvalues.to_markdown())
