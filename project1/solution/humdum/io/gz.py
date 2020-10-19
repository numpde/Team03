# RA, 2020-10-13

import contextlib
import io

@contextlib.contextmanager
def open_maybe_gz(file, *, mode='r'):
    assert mode in ['r', 'rb']

    if isinstance(file, io.IOBase):
        yield file
        return

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
