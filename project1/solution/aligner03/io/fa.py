# RA, 2020-10-09

import typing


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
    """
    try:
        file.seek(0)
    except AttributeError:
        with open(file, mode="r") as file:
            yield from from_fasta(file)
            return

    get_line = (lambda: file.readline().rstrip() or None)

    MARKER = ">"
    sequence = None

    while 1:
        line = get_line()
        if not line:
            assert sequence is not None
            yield sequence
            break
        if line.startswith(MARKER):
            if sequence is not None:
                yield sequence
            sequence = Sequence(line[1:], "")
        else:
            sequence.seq += line
