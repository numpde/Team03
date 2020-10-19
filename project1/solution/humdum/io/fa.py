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
    Yields sequences from FASTA file.
    `file` is a filename, possibly .gz, or io descriptor.

    Assumes there is only one record in the FASTA file.
    """

    from humdum.io import open_maybe_gz

    with open_maybe_gz(file) as file:
        MARKER = ">"

        header = file.readline().strip()
        assert header.startswith(MARKER)
        sequence = Sequence(header[1:], "")
        sequence.seq = "".join(map(str.strip, file.readlines()))

        assert (MARKER not in sequence.seq), "Can only read one record from FASTA."

        yield sequence
