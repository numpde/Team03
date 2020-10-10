# RA, 2020-10-09

import aligner03.io
import collections
import typing
import inclusive
import numpy

_random_state = numpy.random.RandomState(0)


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


def random_kmers(read: aligner03.io.Read, k) -> typing.Iterable[typing.Tuple[int, str, float]]:
    """
    Generate all k-mers of length `k` from the `read` in random order.

    Yields tuples (loc, kmer, qual) where
        `loc` is the location in the read (0-based),
        `kmer` is the kmer,
        `qual` is the average phred score.
    """
    seq = read.seq
    phred = read.phred
    ii = numpy.arange(len(read) - k)
    _random_state.shuffle(ii)
    for i in ii:
        yield (i, seq[i:(i + k)], numpy.mean(phred[i:(i + k)]))


def propose_window(*, read_length: int, read_loc: int, ref_length: int, ref_loc: int) -> typing.Tuple[int, int]:
    """
    Returns a tuple (a, b) such that
        reference_genome[a:b]
    is a piece of genome where to look for the read.

    The window will have roughly twice the size of the read.
    """

    grace_margin = read_length // 2
    a = int(ref_loc - read_loc - grace_margin)
    b = int(a + read_length + (2 * grace_margin))
    a = max(a, 0)
    b = min(b, ref_length)

    return (a, b)
