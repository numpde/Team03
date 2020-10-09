# RA, 2020-10-09

import pysam
import typing

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
