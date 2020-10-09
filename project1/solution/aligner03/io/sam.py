# RA, 2020-10-09

import pysam
import typing
import pandas
import collections


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
            yield read


class FlagFormat:
    @staticmethod
    def has(flag: int, component: int) -> bool:
        return (flag & component) == component

    @staticmethod
    def isMinusStrand(flag):
        return FlagFormat.has(flag, 16)


# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2017/Day1/Session5-alignedReads.html
pandas.Series(name="Flag", dtype=int, data=dict([
    tuple(l.strip().split('\t')) for l in (
        """
        isPaired	1
        isProperPair	2
        isUnmappedQuery	4
        hasUnmappedMate	8
        isMinusStrand	16
        isMateMinusStrand	32
        isFirstMateRead	64
        isSecondMateRead	128
        isSecondaryAlignment	256
        isNotPassingQualityControls	512
        isDuplicate	1024
        """
    ).split('\n') if len(l.strip())
]))
