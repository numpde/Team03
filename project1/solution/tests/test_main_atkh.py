# RA, 2020-10-10

from unittest import TestCase
from pathlib import Path

from humdum.main import AllTheKingsHorses
from humdum.utils import relpath, unlist1

from humdum.io import AlignedSegment

from humdum.io import from_sam_pysam
from pysam import AlignedSegment as pysam_AlignedSegment

from itertools import count


class TestATKH(TestCase):
    def test_on_data_small(self):
        data_root = Path(__file__).parent / "data_for_tests"
        source_path = data_root / "data_small"

        (read_file1, read_file2) = sorted(source_path.glob("*.fq"))
        genome_file = unlist1(source_path.glob("genome*.fa"))

        aligned_segments = AllTheKingsHorses.from_files(fa=genome_file, fq1=read_file1, fq2=read_file2)

        mine: AlignedSegment
        theirs: pysam_AlignedSegment
        for ((mine, theirs), n) in zip(zip(aligned_segments, from_sam_pysam(unlist1(source_path.glob("*.sam")))), count()):
            # See io/sam.py for the explanations
            self.assertEqual(mine.flag.is_minus_strand, bool(theirs.flag & 16))
            self.assertEqual(mine.flag.is_secondary_alignment, bool(theirs.flag & 256))

            # NOTE: pysam returns zero-based indexes
            adjust_pysam_index = (lambda pos: pos + 1)

            cigar_match = (mine.cigar == theirs.cigarstring)
            pos_match = (mine.pos == adjust_pysam_index(theirs.pos))

            if cigar_match and pos_match:
                print(F"Read {mine.qname} looks good.")
            else:
                print(F"Read {mine.qname} does not match.")
                print(F"Mine:  ", mine.cigar, "at", mine.pos)
                print(F"Theirs:", theirs.cigarstring, "at", adjust_pysam_index(theirs.pos))
                print(F"Read:  ", mine.seq)
                #print(F"Neighborhood:  ", atkh.ref_genome[(mine.pos - 10):(mine.pos + 10 + len(mine.seq))])


    def test_header(self):
        AllTheKingsHorses.header()
