# RA, 2020-10-19

from idiva import log

import typing
import contextlib
import io


@contextlib.contextmanager
def seek_then_rewind(fd: io.IOBase, seek=0) -> typing.IO:
    pos = fd.tell()
    if seek is not None:
        fd.seek(seek)
    try:
        yield fd
    finally:
        fd.seek(pos)
