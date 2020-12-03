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


# Adapted from
# https://stackoverflow.com/a/4213255
def checksum_md5(fd: io.TextIOBase) -> str:
    log.info("Computing md5 hash of stream.")
    import hashlib
    with seek_then_rewind(fd):
        md5 = hashlib.md5()
        try:
            for chunk in iter(lambda: fd.read(128 * md5.block_size).encode(), b''):
                md5.update(chunk)
        except KeyboardInterrupt:
            pass
    md5 = md5.hexdigest()
    log.debug(F"Computed md5 = {md5}.")
    return md5
