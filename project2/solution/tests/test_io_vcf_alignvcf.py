# RA, 2020-11-14

import io

from unittest import TestCase
from pathlib import Path

# LARGE DATA EXTRACT

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control_v2.vcf",
    'case': BASE / "case_processed_v2.vcf",
}

# LARGE DATA

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


# TESTS

class TestAlignVCF(TestCase):
    def test_align_largehead(self):
        from idiva.io.vcf import ReadVCF
        from idiva.io.vcf import align

        with PATHS['case'].open() as fd_case:
            with PATHS['ctrl'].open() as fd_ctrl:
                assert isinstance(fd_case, io.TextIOBase)
                assert isinstance(fd_ctrl, io.TextIOBase)
                df = align(case=ReadVCF(fd_case), ctrl=ReadVCF(fd_ctrl))

    def test_align_large(self):
        from idiva.io.vcf import ReadVCF
        from idiva.io.vcf import align
        from idiva.io import open_maybe_gz

        with open_maybe_gz(URLS['case'], mode='r') as fd_case:
            with open_maybe_gz(URLS['ctrl'], mode='r') as fd_ctrl:
                assert isinstance(fd_case, io.TextIOBase)
                assert isinstance(fd_ctrl, io.TextIOBase)
                df = align(case=ReadVCF(fd_case), ctrl=ReadVCF(fd_ctrl))
                print(df)
