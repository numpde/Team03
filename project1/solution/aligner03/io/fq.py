# RA, 2020-09-29

# # FASTQ
# @NAME   DESC KSDSKDJBK
# ATCGTACGTAGCT
# +GHFGG
# hdliufahriuth

import typing

PHRED_OFFSET = 0x21

class Read:
    name: str
    desc: str
    seq: str
    # Note: Ignore the "+" line
    phred: typing.List[int]

    def __init__(self, name: str, desc: str, seq: str, phred: typing.List[int], is_forward=True):
        self.name = name
        self.desc = desc
        self.seq = seq
        self.phred = phred
        self.is_forward = is_forward

    def __str__(self):
        return F"Name: {self.name} ({self.desc}): {self.seq} ({self.phred})"

    @property
    def phred_as_string(self):
        return ''.join([chr(i + PHRED_OFFSET) for i in self.phred])

    @property
    def reverse(self):
        from aligner03.utils.strings import reverse as reverse_complement
        return Read(self.name, self.desc, reverse_complement(self.seq), self.phred[-1::-1], not self.is_forward)


def from_fastq(file) -> typing.Iterable[Read]:
    """
    Yields FASTQ reads from file (filename or descriptor).
    Note: seek(0) is called.
    """

    def phred_as_ints(phred_as_string: str) -> typing.List[int]:
        # https://en.wikipedia.org/wiki/FASTQ_format
        phred = [(ord(c) - PHRED_OFFSET) for c in phred_as_string]
        assert (len(phred) == len(phred_as_string))
        return phred

    try:
        file.seek(0)
    except:
        # Assume `file` needs to be opened first
        with open(file, mode="r") as file:
            yield from from_fastq(file)
            return
    else:
        while 1:
            try:
                SEP = '\t'
                line = (lambda: file.readline().rstrip() or None)
                (at_name, desc) = (line().split(SEP, 1) + [None])[0:2]
                seq = line()
                line()
                phred = line()
            except AttributeError:
                # Assume that line().split failed
                break
            else:
                yield Read(at_name[1:].strip(), desc, seq, phred_as_ints(phred))


def example():
    for read in from_fastq(filename):
        print(read)
