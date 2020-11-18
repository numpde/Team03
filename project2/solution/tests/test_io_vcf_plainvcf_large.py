# RA, 2020-11-12

import typing
import io
from contextlib import ExitStack

from unittest import TestCase
from pathlib import Path
from tcga.utils import download

from idiva.io.gz import open_maybe_gz

download_cache = (Path(__file__).parent.parent.parent / "input/download_cache")
assert download_cache.is_dir()
download = download.to(abs_path=download_cache)

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class TestReadBigFile(TestCase):
    def test_count(self):
        from idiva.io.vcf import ReadVCF, RawDataline
        reference = {'ctrl': 2329288, 'case': 2360972}
        for group in URLS:
            data = download(URLS[group]).now

            with data.open(mode='rb') as fd:
                with open_maybe_gz(fd, mode='r') as fd:
                    vcf = ReadVCF(fd)
                    line: RawDataline
                    nlines = sum(1 for __ in vcf.datalines)
                    # print(F"Group {group} has {nlines} datalines")
                    self.assertEqual(reference[group], nlines)
