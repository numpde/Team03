# RA, 2020-10-09

import typing
from itertools import count


class Sequence:
    desc: str
    seq: str

    def __init__(self, desc, seq):
        self.desc = desc
        self.seq = seq


def from_fasta(file) -> typing.Iterable[Sequence]:
    """
    Yields sequences from FASTA file (filename or descriptor).
    Note: seek(0) is called.

    Assumes there is only one record in the FASTA file.
    """
    try:
        file.seek(0)
    except AttributeError:
        from humdum.io import open_maybe_gz
        with open_maybe_gz(file) as file:
            yield from from_fasta(file)
            return

    MARKER = ">"

    header = file.readline().strip()
    assert header.startswith(MARKER)
    sequence = Sequence(header[1:], "")
    sequence.seq = "".join(map(str.strip, file.readlines()))

    assert MARKER not in sequence.seq

    yield sequence
