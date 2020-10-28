# RA, 2020-10-28

import numpy


def mapq_hist(file):
    """
    Returns a structure containing the fields
        mapq
        counts
    where
        counts[i] is the number of
        transcripts of quality mapq[i]

    Expects a SAM file `file`.

    Note: mapq are neither rounded nor clipped.
    Only reads with |tlen| <= 10000 are considered.
    """

    from humdum.io import from_sam
    from collections import Counter

    mapq_and_counts = numpy.asarray(list(Counter([
        read.mapq
        for read in from_sam(file)
        if (abs(read.tlen) <= 10000)
    ]).items())).T

    class _:
        mapq = mapq_and_counts[0]
        counts = mapq_and_counts[1]

    return _
