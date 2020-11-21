# RA, 2020-11-19

import pandas as pd

from itertools import product
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

from idiva.stat import chi2_test
from scipy.stats import chi2_contingency


class TestChi2(TestCase):
    def test_chi2_oneunit(self):
        (a1, a2) = ([5, 11, 20], [10, 11, 20])
        (c1, c2) = (["A1", "B1", "C1"], ["A2", "B2", "C2"])
        df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])

        (__, ref, __, __) = chi2_contingency([a1, a2])
        cdt = chi2_test(df, (c1, c2)).squeeze()

        self.assertEqual(cdt, ref)

    def test_chi2_na(self):
        (a1, a2) = ([5, 11, 20], [10, 11, pd.NA])
        (c1, c2) = (["A1", "B1", "C1"], ["A2", "B2", "C2"])
        df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])
        cdt = chi2_test(df, (c1, c2)).squeeze()
        self.assertTrue(pd.isna(cdt))

    def test_chi2_large_head(self):
        from idiva.io import ReadVCF, open_maybe_gz
        from idiva.clf.df import v0_df, join

        dfs = {}

        for (k, file) in PATHS_LARGE_HEAD.items():
            with open_maybe_gz(file, mode='r') as fd:
                dfs[k] = v0_df(ReadVCF(fd))

        df = join(case=dfs['case'], ctrl=dfs['ctrl'])

        cols = tuple([F"ALT{n}_{kind}" for n in range(3)] for kind in ['case', 'ctrl'])

        pp = chi2_test(df[cols[0] + cols[1]], cols, add=1)

    def test_chi2_large(self):
        from idiva.io import ReadVCF, open_maybe_gz
        from idiva.clf.df import v0_df, join

        dfs = {}

        for k in URLS_LARGE:
            with download(URLS_LARGE[k]).now.open(mode='rb') as fd:
                with open_maybe_gz(fd, mode='r') as fd:
                    dfs[k] = v0_df(ReadVCF(fd))

        df = join(case=dfs['case'], ctrl=dfs['ctrl'])

        cols = tuple([F"ALT{n}_{kind}" for n in range(3)] for kind in ['case', 'ctrl'])

        p = chi2_test(df[cols[0] + cols[1]], cols, add=1)
