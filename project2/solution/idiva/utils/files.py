# RA, 2020-10-19

import typing
import contextlib


@contextlib.contextmanager
def seek_then_rewind(fd: typing.IO, seek=0) -> typing.IO:
    pos = fd.tell()
    fd.seek(seek)
    try:
        yield fd
    finally:
        fd.seek(pos)
