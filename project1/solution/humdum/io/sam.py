# RA, 2020-10-09

import pysam
import typing
import collections


# https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
# The @HD line should be present, with either the SO tag or the GO tag (but not both) specified
# The @SQ lines should be present if reads have been mapped
# @HD File-level metadata.
#   Optional. If present, there must be
#   only one @HD line and it must be the first line of the file.
# @HD VN:1.6 SO:coordinate
# VN:Format version. Accepted format: /^[0-9]+\.[0-9]+$/
# SO:sorting unknown, unsorted, queryname, coordinate
# For coordinate sort, the major sort key is the RNAME field,
#   with order defined by the order of @SQ lines in the header.
#
# @SQ SN:ref LN:45
# see p.4

# Example:
# @HD	VN:1.4	SO:unsorted
# @SQ	SN:22_5K	LN:4980


# https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
# From the problem description:
# Needs to be present:
#   qname
#   flag (is reversed, is secondary alignment)
#   rname
#   pos
#   cigar
# All other fields can be omitted
# (if the specification allows)
# or be set to their default values
# (0 for integers, * for strings).

# On secondary alignment
# https://www.biostars.org/p/206396/


# THIS IS UGLY
class _:
    is_minus_strand = 16
    is_secondary_alignment = 256

    # # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2017/Day1/Session5-alignedReads.html
    # pandas.Series(name="Flag", dtype=int, data={
    #     "isPaired": 1,
    #     "isProperPair": 2,
    #     "isUnmappedQuery": 4,
    #     "hasUnmappedMate": 8,
    #     "isMinusStrand": 16,
    #     "isMateMinusStrand": 32,
    #     "isFirstMateRead": 64,
    #     "isSecondMateRead": 128,
    #     "isSecondaryAlignment": 256,
    #     "isNotPassingQualityControls": 512,
    #     "isDuplicate": 1024,
    # })

    @classmethod
    def getter(cls, f):
        return (lambda flag: bool(flag.value & getattr(cls, f.__name__)))

    @classmethod
    def setter(cls, f):
        def g(flag, b: bool):
            mask = getattr(cls, f.__name__)
            flag.value = flag.value - (flag.value & mask) + (mask * bool(b))

        return g


class Flag:
    def __init__(self, value: int):
        self.value = value

    # def has(self, component: int) -> bool:
    #     return (self.value & component) == component

    @property
    @_.getter
    def is_minus_strand(self):
        return None

    @is_minus_strand.setter
    @_.setter
    def is_minus_strand(self, v: bool):
        pass

    @property
    @_.getter
    def is_secondary_alignment(self):
        return None

    @is_secondary_alignment.setter
    @_.setter
    def is_secondary_alignment(self, v: bool):
        pass

    def __repr__(self):
        return F"Flag({self.value})"


class AlignedSegment:
    def __init__(self):
        self.qname = "*"
        self.flag = Flag(0)
        self.rname = "*"
        self.pos = 0
        self.mapq = 0
        self.cigar = "*"
        self.rnext = "*"
        self.pnext = 0
        self.tlen = 0
        self.seq = "*"
        self.qual = "*"

    def __str__(self):
        return '\t'.join(map(str, [
            self.qname or "*",
            self.flag.value or 0,
            self.rname or "*",
            self.pos or 0,
            self.mapq or 0,
            self.cigar or "*",
            self.rnext or "*",
            self.pnext or 0,
            self.tlen or 0,
            self.seq or "*",
            self.qual or "*",
        ]))


def from_sam(file) -> typing.Iterable[pysam.AlignedSegment]:
    """
    `file` can be file name or file descriptor.
    Note: seek(0) is called.
    """

    try:
        file.seek(0)
    except AttributeError:
        with open(file, mode='r') as file:
            yield from from_sam(file)
            return

    # https://pysam.readthedocs.io/en/latest/api.html
    with pysam.AlignmentFile(file) as af:
        for read in af.fetch():
            read: pysam.AlignedSegment
            yield read
