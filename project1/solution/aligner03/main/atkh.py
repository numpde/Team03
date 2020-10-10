# RA, 2020-10-10

"""
The decision nexus.
"""

from aligner03.index import FmIndex as GenomeIndex
from aligner03.align import SmithWaterman as SequenceAligner
from aligner03.align import Alignment

from aligner03.io import assert_order_consistency, from_fastq, Read
from aligner03.io import AlignedSegment
from aligner03.map import random_kmers, propose_window
from aligner03.utils import at_most_n, first

import typing
import numpy


class UnmappedReadpair(Exception):
    pass


class AllTheKingsHorses:
    """
    Attempts to undo shotgunning.
    """

    class _:
        kmers_per_read = 5
        seed_kmer_size = 26

    def __init__(self, genome_index: GenomeIndex, sequence_aligner: SequenceAligner, ref_genome=None):
        self.index = genome_index
        self.align = sequence_aligner
        self.ref_genome = ref_genome or str(genome_index)

    def map_one(self, read, decide=True) -> typing.Dict[Read, typing.List]:
        """
        The primary purpose of this is to find whether
        the read aligns forward or backward to the reference.

        'Backward' = reverse complement.
        """

        proposals = {
            r: [
                (loc_in_read, None, qual, loc_in_ref)
                for (loc_in_read, kmer, qual) in random_kmers(r, k=self._.seed_kmer_size, maxn=self._.kmers_per_read)
                for loc_in_ref in list(self.index.query(kmer))
            ]
            for r in [read, read.reversed]
        }

        if decide:
            # Keep the read with more matches
            return dict([max(proposals.items(), key=(lambda p: len(p[1])))])
        else:
            return proposals

    def select_option(self, options):
        options: typing.List[typing.Tuple[int, str, float, int]]
        in_ref = [j for (_, _, _, j) in options]
        j = int(numpy.percentile(in_ref, 50, interpolation='nearest'))
        return max(options, key=(lambda x: (x[3] == j)))

    def map_pair(self, read1, read2):
        (read1, options1) = self.map_one(read1).popitem()
        (read2, options2) = self.map_one(read2).popitem()

        if (read1.is_forward == read2.is_forward):
            raise UnmappedReadpair

        if (not options1) or (not options2):
            raise UnmappedReadpair

        ref_length = len(self.ref_genome)

        (i1, _, _, j1) = self.select_option(options1)
        (i2, _, _, j2) = self.select_option(options2)

        w1 = propose_window(read_length=len(read1), read_loc=i1, ref_length=ref_length, ref_loc=j1)
        w2 = propose_window(read_length=len(read2), read_loc=i2, ref_length=ref_length, ref_loc=j2)

        w1_segment = self.ref_genome[w1[0]:w1[1]]
        w2_segment = self.ref_genome[w2[0]:w2[1]]

        alignment1: Alignment
        alignment2: Alignment
        alignment1 = first(self.align(ref=w1_segment, query=read1.seq))
        alignment2 = first(self.align(ref=w2_segment, query=read2.seq))

        alignment1._start_pos += w1[0]
        alignment2._start_pos += w2[0]


        for (read, alignment) in zip([read1, read2], [alignment1, alignment2]):
            seg = AlignedSegment()
            seg.qname = read.preprocessed.name
            seg.flag.is_minus_strand = read.is_forward # Correct?
            seg.cigar = alignment.cigar
            print(alignment._start_pair)
            # print(seg.pos)
            # print(seg)


        # Now we have one read forward and one backward
        # Expect that forward comes before backward

    def map_paired(self, file1, file2):
        assert_order_consistency(file1, file2)

        unmapped_pairs = 0

        for (read1, read2) in zip(from_fastq(file1), from_fastq(file2)):
            try:
                yield self.map_pair(read1, read2)
            except UnmappedReadpair:
                unmapped_pairs += 1
                print("Unmapped reads:", unmapped_pairs)

    def print_header(self):
        # TODO
        raise NotImplementedError
