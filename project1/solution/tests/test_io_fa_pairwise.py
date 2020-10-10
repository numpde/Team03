# RA, 2020-10-10

"""
Check that the reads are in the same order
in FASTQ files.
"""

from pathlib import Path
from unittest import TestCase
from aligner03.io import assert_order_consistency


class TestPairwise(TestCase):
    def test_is_pairwise_data_small(self):
        (file1, file2) = sorted((Path(__file__).parent / "data_for_tests/data_small/").glob("*.fq"))
        assert_order_consistency(file1, file2)

    def test_is_pairwise_data_large(self):
        (file1, file2) = sorted((Path(__file__).parent / "data_for_tests/data/").glob("*.fq"))
        assert_order_consistency(file1, file2)
