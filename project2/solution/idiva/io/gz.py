# RA, 2020-10-13

import contextlib
import io
import typing

@contextlib.contextmanager
def open_maybe_gz(file, *, mode='r') -> typing.IO:
    """
    Open `file` for reading that could be a
     - file descriptor
     - path to file
     - path to gzipped file

    `mode` is either 'r' or 'rb', and has to be specified.

    Usage:
        with open_maybe_gz(path_to_file, mode='r') as fd:
            print(fd.read())
    """

    assert mode in ['r', 'rb']

    if isinstance(file, io.IOBase):
        yield file
        return

    from pathlib import Path
    assert Path(file).is_file()

    file = str(file)

    if file.endswith(".gz"):
        import gzip
        with gzip.open(file, mode='rb') as fd:
            if (mode == 'r'):
                yield io.TextIOWrapper(fd)
            elif (mode == 'rb'):
                yield fd
    else:
        with open(file, mode=mode) as fd:
            yield fd
