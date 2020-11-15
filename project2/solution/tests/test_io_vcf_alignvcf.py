# RA, 2020-11-14

import pandas as pd
import numpy as np
import typing
import io

from unittest import TestCase
from pathlib import Path
from idiva.utils import unlist1, at_most_n, first
from tcga.utils import download

# LARGE DATA EXTRACT

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control.vcf",
    'case': BASE / "case_processed.vcf",
}

# LARGE DATA

download_cache = (Path(__file__).parent.parent.parent / "input/download_cache")
assert download_cache.is_dir()
download = download.to(abs_path=download_cache)

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control.vcf",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed.vcf",
}


# TESTS

class TestAlignVCF(TestCase):
    def test_align_largehead(self):
        from idiva.io.vcf import ReadVCF
        from idiva.io.vcf import align

        with PATHS['case'].open() as fd_case:
            with PATHS['ctrl'].open() as fd_ctrl:
                df = align(case=ReadVCF(fd_case), ctrl=ReadVCF(fd_ctrl))
                print(df)

    def test_align_large(self):
        from idiva.io.vcf import ReadVCF
        from idiva.io.vcf import align

        data = {k: download(URLS[k]).now for k in URLS}
        with data['case'].open() as fd_case:
            with data['ctrl'].open() as fd_ctrl:
                df = align(case=ReadVCF(fd_case), ctrl=ReadVCF(fd_ctrl))
                print(df)
