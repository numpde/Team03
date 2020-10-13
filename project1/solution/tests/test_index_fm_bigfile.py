# RA, LB

from unittest import TestCase
from pathlib import Path

from humdum.io import open_maybe_gz
from humdum.io import from_fasta
from humdum.utils import unlist1, first

data_root = Path(__file__).parent / "data_for_tests/data"
genome_file = unlist1(data_root.glob("*.fa.gz"))


class TestFm(TestCase):
    def test_open_and_read(self):
        with open_maybe_gz(genome_file) as fd:
            (fd.readline())
            (fd.readline())
            (fd.readline())

    def test_index(self):
        genome = first(from_fasta(genome_file))
        print(len(genome))
        raise NotImplementedError
