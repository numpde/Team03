# RA, 2020-10-09

import aligner03.io
import collections
import typing
import inclusive
import numpy


def all_kmers_by_score(read: aligner03.io.Read, k: int) -> typing.Dict[float, typing.List[tuple]]:
    """
    Get all kmers arranged by phred score.
    """
    by_score = collections.defaultdict(list)
    for i in inclusive.range[0, len(read.seq) - k]:
        ii = slice(i, i + k)
        by_score[
            numpy.average(read.phred[ii])
        ].append(
            (i, read.seq[ii])
        )
    return dict(by_score)
