# RA, 2020-11-14

import typing
import io

from unittest import TestCase
from pathlib import Path
from idiva.utils import unlist1, at_most_n, first
from tcga.utils import download

download_cache = (Path(__file__).parent.parent.parent / "input/download_cache")
assert download_cache.is_dir()
download = download.to(abs_path=download_cache)

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control.vcf",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed.vcf",
}

# UNFINISHED

class TestAlignVCF(TestCase):
    def test_align(self):
        from idiva.io.vcf import ReadVCF, RawDataline
        from idiva.io.vcf import AlignVCF

        data = {k: download(URLS[k]).now for k in URLS}
        with data['ctrl'].open() as fd1, data['case'].open() as fd2:
            vcf1 = ReadVCF(fd1)
            vcf2 = ReadVCF(fd2)
            for (n, (d1, d2)) in enumerate(AlignVCF(vcf1, vcf2), start=1):
                print(n)
                print(d1)
                print(d2)

                self.assertEqual(d1.pos, d2.pos)
                self.assertEqual(d1.id, d2.id)


