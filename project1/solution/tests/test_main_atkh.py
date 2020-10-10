# RA, 2020-10-10

from unittest import TestCase
from pathlib import Path
from aligner03.main import AllTheKingsHorses
from aligner03.utils import relpath, unlist1

class TestATKH(TestCase):
    def test_on_data_small(self):
        data_root = Path(__file__).parent / "data_for_tests"
        source_path = data_root / "data_small"
        (file1, file2) = sorted(source_path.glob("*.fq"))


        genome_file_fasta = unlist1(source_path.glob("genome*.fa"))
        
        AllTheKingsHorses()