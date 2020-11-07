# RA, 2020-11-05

import typing


class Oneliner:
    def __init__(self, fd: typing.TextIO):
        self.fd = fd

    @property
    def last(self):
        assert hasattr(self, "line"), "You need to start iterating over me first."
        return self.line

    def __next__(self):
        while 1:
            line = self.fd.readline()
            if line:
                self.line = line.strip("\n")
                if self.line:
                    return self.line
            else:
                raise StopIteration

    def __iter__(self):
        return self
