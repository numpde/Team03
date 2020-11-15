# RA, 2020-11-05

import typing


class Oneliner:
    """
    Read a TextIO line by line but skip empty ones.
    Trims the lines using str.strip.

    for line in Oneliner(fd):
        print(line)

    The attribute `last` keeps the last extracted line.
    """

    def __init__(self, fd: typing.TextIO, buffered=True):
        self.fd = fd
        self.buffered = buffered
        self._line_stack = []

    @property
    def last(self):
        assert hasattr(self, 'line'), "You need to start iterating over me first."
        return self.line

    def __next__(self):
        while 1:
            if not self._line_stack:
                if self.buffered:
                    self._line_stack = list(map(str.strip, reversed(self.fd.readlines(1024))))
                else:
                    line = self.fd.readline()
                    if line:
                        self._line_stack = [line]
            if not self._line_stack:
                raise StopIteration
            self.line = self._line_stack.pop()
            # Skip empty lines
            if self.line:
                return self.line

    def __iter__(self):
        return self
