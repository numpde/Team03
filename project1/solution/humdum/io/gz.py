# RA, 2020-10-13

import contextlib
import io

@contextlib.contextmanager
def open_maybe_gz(file):
    file = str(file)
    if file.endswith(".gz"):
        import gzip
        with gzip.open(file, mode='r') as fd:
            yield io.TextIOWrapper(fd)
    else:
        with open(file, mode='r') as fd:
            yield fd
