# RA, 2020-10-13

from pathlib import Path
from unittest import TestCase
from itertools import count

from humdum.utils import unlist1, at_most_n

data_root = Path(__file__).parent / "data_for_tests"
source_path = data_root / "data_small"
file = unlist1(sorted(source_path.glob("*.sam")))


class TestIoSam(TestCase):
    def test_import(self):
        from humdum.io import from_sam
        pass

    def test_no_import(self):
        """
        Test whether from_sam_pysam is gone.
        """
        with self.assertRaises(ImportError):
            from humdum.io import from_sam_pysam
            pass

    def test_from_sam(self):
        from humdum.io import from_sam

        alignment = list(from_sam(file)).pop()

        last_read = "CACCATCCAGAACAGTGCCTCTTGCAGAGTCTCCTTGGGAAACTTACCAAGTCTGATGGTAGCAGGGGCATGGGACCATCCTAACTGGGAAGACAAAAAGGCTGAGACCTTCCCAGAGTCACCTT"
        self.assertEqual(alignment.seq, last_read)


    def test_against_pysam(self):
        from humdum.io import from_sam, AlignedSegment
        from humdum.io import from_sam_pysam
        import pysam

        for (have, want, n) in zip(from_sam(file), from_sam_pysam(file), count()):
            self.assertIsInstance(have, AlignedSegment)
            self.assertIsInstance(want, pysam.AlignedSegment)
            self.assertEqual(have.cigar, want.cigarstring)

        self.assertEqual(n, 1169)
