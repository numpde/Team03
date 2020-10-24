# RA, 2020-10-24

import numpy


def tlen_hist(file):
    """
    Return a structure containing the fields
        length
        counts
    where
        counts[i] is the number of
        transcripts of length length[i]

    Expectes a SAM file `file`.

    Only counts nonnegative tlen.
    """

    from humdum.io import from_sam
    from collections import Counter

    tlens_counts = numpy.asarray(list(Counter([
        read.tlen
        for read in from_sam(file)
        if (read.tlen >= 0)
    ]).items())).T

    class _:
        length = tlens_counts[0]
        counts = tlens_counts[1]

    return _
