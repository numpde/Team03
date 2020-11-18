# RA, 2020-11-14

"""
Check the assumptions on VCF files from idiva.io.ass .
"""

import typing
import io

from unittest import TestCase
from pathlib import Path

from idiva.io import open_maybe_gz
from idiva.utils import unlist1, at_most_n, first
from tcga.utils import download

download_cache = (Path(__file__).parent.parent.parent / "input/download_cache")
assert download_cache.is_dir()
download = download.to(abs_path=download_cache)

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestAssumptions(TestCase):
    def test_check_all(self):
        for group in ['case', 'ctrl']:
            data = download(URLS[group]).now
            with data.open(mode='rb') as fd:
                with open_maybe_gz(fd) as fd:
                    from idiva.io.ass import check_all
                    check_all(fd)
