# RA, 2020-09-28

try:
    import pysam

    Read = pysam.libcalignedsegment.AlignedSegment
except:
    print("Need the package pysam")
    raise

import numpy as np

from pathlib import Path
from typing import Iterable


# On alignment quality

# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5964631/
# Three primary metrics are used to evaluate sequence quality in NGS data:
# depth of coverage (how many sequence reads are present at a given position),
# base quality (have the correct bases been called in sequence reads) and
# mapping quality (have the reads been mapped to the correct position in the genome).

# https://learn.gencore.bio.nyu.edu/alignment/


def read_sam(file) -> Iterable[Read]:
    """
    `file` can be file name or file descriptor
    """

    if (type(file) is str):
        with open(file, mode='r') as file:
            yield from read_sam(file)
        return

    # https://pysam.readthedocs.io/en/latest/api.html
    with pysam.AlignmentFile(file) as af:
        for read in af.fetch():
            yield read


def coverage_pbp(file, reference_length=None):
    """
    Reads the SAM file and computes the per-based coverage
    from the aligned blocks.
    In particular: if a position in the reference is consistently
    a 'delete' compared to the matched reads, it will be counted as zero.

    The `reference_length` can be inferred from the mapped reads
    but there can be an undetected residual at the 'right' end
    if it is not specified.
    """
    zeros = (lambda n: np.zeros(n, dtype=int))
    counts = zeros(reference_length or 0)
    for read in read_sam(file):
        for (a, b) in read.get_blocks():
            # A block is a no-gap alignment
            # Note: Blocks are 1- based
            assert (a < b)
            if (b > len(counts)):
                counts = np.concatenate([counts, zeros(b - len(counts))])
            counts[(a - 1):b] += 1

    return counts


def test_read_sam(file):
    with open(file, mode='rb') as fd_sam:
        for read in read_sam(fd_sam):
            read: Read

            # https://en.wikipedia.org/wiki/SAM_(file_format)
            # Col   Field	Type	Brief description
            #   1   QNAME	String  Query template NAME
            #   2   FLAG	Int	    bitwise FLAG
            #   3   RNAME	String  References sequence NAME
            #   4   POS     Int     1- based leftmost mapping POSition
            #   5   MAPQ	Int     MAPping Quality
            #   6   CIGAR	String	CIGAR String
            #   7   RNEXT	String	Ref. name of the mate/next read
            #   8   PNEXT	Int     Position of the mate/next read
            #   9   TLEN	Int     observed Template LENgth
            #  10	SEQ     String  segment SEQuence
            #  11	QUAL	String  ASCII of Phred-scaled base QUALity+33

            print("QNAME ", read.query_name)
            print("FLAG  ", read.flag)
            print("RNAME ", read.reference_id)
            print("POS   ", read.reference_start)
            print("MAPQ  ", read.mapping_quality)
            print("CIGAR ", read.cigartuples)
            print("RNEXT ", read.next_reference_id)
            print("PNEXT ", read.next_reference_start)
            print("TLEN  ", read.template_length)
            print("SEQ   ", read.query_sequence)
            print("QUAL  ", read.query_qualities)

            # print(read.get_aligned_pairs())
            print(read.get_blocks())


def test_coverage_pbp(file):
    counts = coverage_pbp(file)
    print(*counts)


if (__name__ == '__main__'):
    in_file = list((Path(__file__).parent.parent.parent.parent / "input/data_small/").glob("*.sam")).pop()
    # test_read_sam(in_file)
    test_coverage_pbp(in_file)
