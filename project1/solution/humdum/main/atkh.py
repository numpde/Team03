# RA, 2020-10-10

"""
The decision nexus.
"""

import sys

from humdum.io import from_fasta
from humdum.io import AlignedSegment
from humdum.io import Sequence as FastaGenome
from humdum.io import assert_order_consistency, from_fastq, Read

try:
    from humdum.index import FmIndex as GenomeIndex
    from humdum.qc import healthcheck_index

    healthcheck_index(GenomeIndex)
except Exception as ex:
    print(F"Warning: FmIndex import failed in {__file__} ({ex}).", file=sys.stderr)
else:
    from humdum.index.naive import NaiveIndex as GenomeIndex

from humdum.align import Alignment
from humdum.align import SmithWaterman as SequenceAligner

from humdum.map import random_kmers, propose_window
from humdum.utils import unlist1, first

import typing
import numpy


class UnmappedReadpair(Exception):
    def __init__(self, info):
        self.info = info


class AllTheKingsHorses:
    """
    Attempts to undo shotgunning.
    """

    class Settings:
        kmers_per_read = 5
        seed_kmer_size = 26

    def __init__(self, genome_index: GenomeIndex, sequence_aligner: SequenceAligner, ref_genome: FastaGenome):
        assert isinstance(ref_genome, FastaGenome)

        self.index = genome_index
        self.align = sequence_aligner
        self.ref_genome = ref_genome

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
                for (loc_in_read, kmer, qual) in
                random_kmers(r, k=self.Settings.seed_kmer_size, maxn=self.Settings.kmers_per_read)
                for loc_in_ref in list(self.index.query(kmer))
            ]
            for r in [read, read.reversed]
        }

        if decide:
            # Keep the read with more matches
            return dict([max(proposals.items(), key=(lambda p: len(p[1])))])
        else:
            return proposals

    @staticmethod
    def select_option(options):
        options: typing.List[typing.Tuple[int, str, float, int]]
        in_ref = [j for (_, _, _, j) in options]
        j = int(numpy.percentile(in_ref, 50, interpolation='nearest'))
        return max(options, key=(lambda x: (x[3] == j)))

    def map_pair(self, read1, read2) -> typing.Iterable[AlignedSegment]:
        (read1, options1) = self.map_one(read1).popitem()
        (read2, options2) = self.map_one(read2).popitem()

        if (read1.is_forward == read2.is_forward):
            raise UnmappedReadpair({'reads': [read1, read2], 'reason': "Directionality"})

        if (not options1) or (not options2):
            raise UnmappedReadpair({'reads': [read1, read2], 'reason': "No mapping options"})

        ref_length = len(self.ref_genome.seq)

        read2seg = {}

        # Transcript beginning & end, inclusive range [ta, tb], 0-based
        ta = None
        tb = None

        for (read, options) in zip([read1, read2], [options1, options2]):
            (i, _, _, j) = self.select_option(options)

            w = propose_window(read_length=len(read), read_loc=i, ref_length=ref_length, ref_loc=j)

            w_segment = self.ref_genome.seq[w[0]:w[1]]

            alignment: Alignment
            alignment = first(self.align(ref=w_segment, query=read.seq, alignment_type='semi-local'))

            loc_in_ref = (alignment.loc_in_ref + w[0])

            # Infer transcript extent [ta, tb]
            ta = min(loc_in_ref, loc_in_ref + alignment.tlen - 1, ta or ref_length)
            tb = max(loc_in_ref, loc_in_ref + alignment.tlen - 1, tb or 1)

            seg = AlignedSegment()
            seg.qname = read.preprocessed.name
            # Need to set two flags:
            # is_reversed, is_secondary_alignment
            seg.flag.is_minus_strand = not read.is_forward
            seg.flag.is_secondary_alignment = False
            seg.mapq = alignment.score  # TODO: is this OK
            seg.cigar = alignment.cigar
            seg.pos = loc_in_ref + 1
            seg.seq = read.seq
            seg.qual = read.phred_as_string

            read2seg[read] = seg

        # Get position of mate
        read2seg[read1].pnext = read2seg[read2].pos
        read2seg[read2].pnext = read2seg[read1].pos

        # Get fragment length (https://www.biostars.org/p/356811/)
        tlen = tb - ta + 1  # Absolute length
        for read in [read1, read2]:
            read2seg[read].tlen = -tlen if read2seg[read].flag.is_minus_strand else tlen

        # Yield in the correct order
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
            except UnmappedReadpair as ex:
                self.unmapped_readpairs += 1

    def headers(self) -> typing.Iterable[str]:
        # https://samtools.github.io/hts-specs/SAMv1.pdf
        yield '\t'.join(["@HD", "VN:1.1", "SO:unsorted"])
        yield '\t'.join(["@SQ", "SN:{}".format(self.ref_genome.desc), "LN:{}".format(len(self.ref_genome.seq))])
        yield '\t'.join(["@PG", "ID:humdum"])

    @classmethod
    def from_files(cls, *, fa, fq1, fq2):
        """
        Reference genome file `fa`.
        FASTQ files `fq1` and `fq2`.

        Creates an instance of AllTheKingsHorses and
        yields from its map_paired(...) member function.
        """

        ref_genome = unlist1(from_fasta(fa))

        index = GenomeIndex.read_or_make(path_to_genome=fa)

        aligner = SequenceAligner()

        atkh = AllTheKingsHorses(genome_index=index, sequence_aligner=aligner, ref_genome=ref_genome)

        class _:
            headers = atkh.headers()
            alignments = atkh.map_paired(fq1, fq2)

        return _
