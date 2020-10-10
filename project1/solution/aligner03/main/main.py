# RA, 2020-10-10

"""
The decision nexus.
"""

from aligner03.index import FmIndex as GenomeIndex
from aligner03.align import SmithWaterman as SequenceAligner

from aligner03.io import assert_order_consistency, from_fastq
from aligner03.map import random_kmers, propose_window
from aligner03.utils import at_most_n

class AllTheKingsHorses:
    """
    Attempts to undo shotgunning.
    """

    def __init__(self, genome_index, sequence_aligner, seed_size):
        self.index = genome_index
        self.align = sequence_aligner
        self.k = seed_size

    def map_paired(self, reads_file1, reads_file2):
        assert_order_consistency(reads_file1, reads_file2)

        for (read1, read2) in zip(from_fastq(reads_file1), from_fastq(reads_file2)):
            yield from self.map_pair(read1, read2)

    def map_pair(self, read1, read2):
        self.index: GenomeIndex
        self.align: SequenceAligner

        for (i, kmer, qual) in random_kmers(read1, self.k):
            alignments = list(self.index.query(kmer))
            if alignments:
                print(alignments)
                propose_window(read_length=len(read1), read_loc=i, ref_length=len(self.index), ref_loc=)

        self.index
        yield None

