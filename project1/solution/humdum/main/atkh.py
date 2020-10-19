# RA, 2020-10-10

"""
The decision nexus.
"""

from humdum.io import from_fasta
from humdum.io import AlignedSegment
from humdum.io import assert_order_consistency, from_fastq, Read

from humdum.index import FmIndex as GenomeIndex

from humdum.align import Alignment
from humdum.align import SmithWaterman as SequenceAligner

from humdum.map import random_kmers, propose_window
from humdum.utils import unlist1, at_most_n, first

import typing
import numpy


class UnmappedReadpair(Exception):
    pass


class AllTheKingsHorses:
    """
    Attempts to undo shotgunning.
    """

    class _:
        # Settings
        kmers_per_read = 5
        seed_kmer_size = 26

    def __init__(self, genome_index: GenomeIndex, sequence_aligner: SequenceAligner, ref_genome=None):
        self.index = genome_index
        self.align = sequence_aligner
        self.ref_genome = ref_genome or str(genome_index)
        self.unmapped_readpairs = None

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

    def map_pair(self, read1, read2) -> typing.Iterable[AlignedSegment]:
        (read1, options1) = self.map_one(read1).popitem()
        (read2, options2) = self.map_one(read2).popitem()

        if (read1.is_forward == read2.is_forward):
            raise UnmappedReadpair

        if (not options1) or (not options2):
            raise UnmappedReadpair

        ref_length = len(self.ref_genome)

        read2seg = {}

        for (read, options) in zip([read1, read2], [options1, options2]):
            (i, _, _, j) = self.select_option(options)

            w = propose_window(read_length=len(read), read_loc=i, ref_length=ref_length, ref_loc=j)

            w_segment = self.ref_genome[w[0]:w[1]]

            alignment: Alignment
            alignment = first(self.align(ref=w_segment, query=read.seq, alignment_type='semi-local'))

            loc_in_ref = (alignment.loc_in_ref + w[0])

            seg = AlignedSegment()
            seg.qname = read.preprocessed.name
            # Need to set two flags:
            # is_reversed, is_secondary_alignment
            seg.flag.is_minus_strand = not read.is_forward
            seg.flag.is_secondary_alignment = False
            seg.cigar = alignment.cigar
            seg.pos = loc_in_ref + 1
            seg.seq = read.seq
            seg.qual = read.phred_as_string

            read2seg[read] = seg

        # Get position of mate
        read2seg[read1].pnext = read2seg[read2].pos
        read2seg[read2].pnext = read2seg[read1].pos

        for read in [read1, read2]:
            yield read2seg[read]

    def map_paired(self, file1, file2) -> typing.Iterable[AlignedSegment]:
        """
        Yield aligned segments.

        Check the member variable `unmapped_readpairs` for the number of discarded pairs.
        """

        assert_order_consistency(file1, file2)

        self.unmapped_readpairs = 0

        for (read1, read2) in zip(from_fastq(file1), from_fastq(file2)):
            try:
                yield from self.map_pair(read1, read2)
            except UnmappedReadpair:
                self.unmapped_readpairs += 1

    @classmethod
    def header(self):
        # TODO
        raise NotImplementedError

    @classmethod
    def from_files(cls, *, fa, fq1, fq2):
        """
        Reference genome file `fa`.
        FASTQ files `fq1` and `fq2`.

        Creates an instance of AllTheKingsHorses and
        yields from its map_paired(...) member function.
        """

        ref_genome = unlist1(from_fasta(fa)).seq

        index = GenomeIndex(ref_genome)
        aligner = SequenceAligner()

        atkh = AllTheKingsHorses(genome_index=index, sequence_aligner=aligner, ref_genome=ref_genome)

        for alignment in atkh.map_paired(fq1, fq2):
            yield alignment
