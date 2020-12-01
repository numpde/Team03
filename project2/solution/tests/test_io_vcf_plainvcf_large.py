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
        # ref_len_v1 = {'ctrl': 2329288, 'case': 2360972}
        ref_len_v2 = {'ctrl': 2227080, 'case': 2258797}
        for group in URLS:
            with download(URLS[group]).now.open(mode='rb') as fd:
                with open_maybe_gz(fd, mode='r') as fd:
                    assert isinstance(fd, io.TextIOBase)
                    nlines = sum(1 for __ in ReadVCF(fd))
                    # print(F"Group {group} has {nlines} datalines")
                    self.assertEqual(nlines, ref_len_v2[group])
