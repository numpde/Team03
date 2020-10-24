# RA, 2020-10-23

from unittest import TestCase
from pathlib import Path

from humdum.main import AllTheKingsHorses
from humdum.utils import relpath, unlist1, at_most_n

from humdum.io import AlignedSegment
from humdum.io import from_sam

from itertools import count

data_root = Path(__file__).parent / "data_for_tests"
source_path = data_root / "data"
genome_file = unlist1(source_path.glob("genome*.fa.gz"))


class TestATKH(TestCase):
    def test_on_data_large_5xCov(self):
        (read_file1, read_file2) = sorted(source_path.glob("*5xCov*.fq*"))

        sam = AllTheKingsHorses.from_files(fa=genome_file, fq1=read_file1, fq2=read_file2)

        for alignment in at_most_n(sam.alignments, 50):
            print(alignment)

