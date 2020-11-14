# RA, 2020-10-19

import typing
import contextlib


@contextlib.contextmanager
def seek0(fd: typing.IO):
    pos = fd.tell()
    fd.seek(0)
    yield fd
    fd.seek(pos)
