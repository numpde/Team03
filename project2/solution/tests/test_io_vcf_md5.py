# RA, 2020-12-02

import io
from unittest import TestCase
from pathlib import Path

BASE = (Path(__file__).parent) / "data_for_tests"

PATHS_LARGE_HEAD = {
    'ctrl': BASE / "large_head/control_v2.vcf",
    'case': BASE / "large_head/case_processed_v2.vcf",
}

PATHS_LARGE_FULL = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class Test_ReadVCF_MD5(TestCase):
    def test_md5_head(self):
        from idiva.io import ReadVCF, open_maybe_gz
        from idiva.utils import seek_then_rewind

        for (k, file) in PATHS_LARGE_HEAD.items():
            with open_maybe_gz(file, mode='r') as fd:
                assert isinstance(fd, io.TextIOBase)
                with seek_then_rewind(fd, seek=0):
                    import hashlib
                    reference = hashlib.md5(fd.read().encode()).hexdigest()

                candidate = ReadVCF(fd).md5

                self.assertEqual(candidate, reference)

    def test_md5_full(self):
        from idiva.io import ReadVCF, open_maybe_gz

        for (k, file) in PATHS_LARGE_FULL.items():
            with open_maybe_gz(file, mode='r') as fd:
                assert isinstance(fd, io.TextIOBase)
                ReadVCF(fd).md5

