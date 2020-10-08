# RA, 2020-09-29

# # FASTQ
# @NAME   DESC KSDSKDJBK
# ATCGTACGTAGCT
# +GHFGG
# hdliufahriuth

import typing


class Read:
    name: str
    desc: str
    seq: str
    # Note: Ignore the "+" line
    phred: typing.List[int]

    def __init__(self, name, desc, seq, phred):
        self.name = name
        self.desc = desc
        self.seq = seq
        self.phred = phred

    def __str__(self):
        return F"Name: {self.name} ({self.desc}): {self.seq} ({self.phred})"


def from_fastq(file) -> typing.Iterable[Read]:
    """
    Yields FASTQ reads from file (filename or descriptor).
    """

    def phred_as_ints(phred: str) -> typing.List[int]:
        # TODO
        phred = [None] * len(phred)
        return phred

    if type(file) is str:
        with open(file, mode="r") as file:
            yield from from_fastq(file)
        return

    while 1:
        try:
            SEP = '\t'
            line = (lambda: file.readline().rstrip())
            (at_name, desc) = ((line()).split(SEP, 1) + (None, ))[0:2]
            seq = line()
            line()
            phred = line()

            yield Read(at_name[1:].strip(), desc, seq, phred_as_ints(phred))
        except:  # todo: appropriate exception
            break


def test_reads():
    for read in from_fastq(filename):
        print(read)
