# RA, 2020-10-10

from unittest import TestCase
from pathlib import Path

from humdum.io import from_fasta
from humdum.main import AllTheKingsHorses
from humdum.utils import relpath, unlist1
from humdum.index import FmIndex as GenomeIndex
from humdum.align import SmithWaterman as Aligner
# from humdum.index import NaiveIndex as GenomeIndex

from itertools import count


class TestATKH(TestCase):
    def test_on_data_small(self):
        data_root = Path(__file__).parent / "data_for_tests"
        source_path = data_root / "data_small"
        (file1, file2) = sorted(source_path.glob("*.fq"))

        genome_file = unlist1(source_path.glob("genome*.fa"))
        ref_genome = unlist1(from_fasta(genome_file)).seq
        index = GenomeIndex(ref_genome)

        aligner = Aligner()

        atkh = AllTheKingsHorses(genome_index=index, sequence_aligner=aligner, ref_genome=ref_genome)

        from humdum.io import AlignedSegment
        from pysam import AlignedSegment as pysam_AlignedSegment
        from humdum.io import from_sam

        from_aether = atkh.map_paired(file1, file2)

        mine: AlignedSegment
        theirs: pysam_AlignedSegment
        for ((mine, theirs), n) in zip(zip(from_aether, from_sam(unlist1(source_path.glob("*.sam")))), count()):
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
                print(F"Neighborhood:  ", ref_genome[(mine.pos - 10):(mine.pos + 10 + len(mine.seq))])
