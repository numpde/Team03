# RA, 2020-10-10

from unittest import TestCase
from pathlib import Path

from aligner03.io import from_fasta
from aligner03.main import AllTheKingsHorses
from aligner03.utils import relpath, unlist1
from aligner03.index import FmIndex as GenomeIndex
from aligner03.align import SmithWaterman as Aligner
# from aligner03.index import NaiveIndex as GenomeIndex


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


        from aligner03.io import AlignedSegment
        from pysam import AlignedSegment as pysam_AlignedSegment
        from aligner03.io import from_sam

        from_aether = atkh.map_paired(file1, file2)

        mine: AlignedSegment
        theirs: pysam_AlignedSegment
        for (mine, theirs) in zip(from_aether, from_sam(unlist1(source_path.glob("*.sam")))):
            self.assertEqual(mine.flag.is_minus_strand, bool(theirs.flag & 16))
            self.assertEqual(mine.flag.is_secondary_alignment, bool(theirs.flag & 256))

            cigar_match = (mine.cigar == theirs.cigarstring)

            if not cigar_match:
                print(F"Read {mine.qname}, has non-matching cigar.")
                print(F"Mine:  ", mine.cigar)
                print(F"Theirs:", theirs.cigarstring)
                print(F"Read:  ", mine.seq)
                print(F"Neighborhood:  ", ref_genome[(mine.pos - 10):(mine.pos + 10 + len(mine.seq))])


