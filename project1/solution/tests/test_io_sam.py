# RA, 2020-10-13

from pathlib import Path
from unittest import TestCase

from humdum.utils import unlist1, at_most_n

data_root = Path(__file__).parent / "data_for_tests"
source_path = data_root / "data_small"


class TestIoSam(TestCase):
    def test_import(self):
        from humdum.io import from_sam
        pass

    def test_no_import(self):
        with self.assertRaises(ImportError):
            from humdum.io import from_sam_pysam
            pass

    def test_from_sam(self):
        from humdum.io import from_sam

        file = unlist1(sorted(source_path.glob("*.sam")))
        sam = from_sam(file)
        for alignment in at_most_n(sam, 10):
            alignment
