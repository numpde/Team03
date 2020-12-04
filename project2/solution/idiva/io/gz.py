# RA, 2020-10-13

import contextlib
import io
import typing
import re

from idiva.utils import seek_then_rewind


@contextlib.contextmanager
def open_maybe_gz(file, *, mode='r') -> typing.Union[io.TextIOBase, io.BytesIO]:
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

    if isinstance(file, io.TextIOBase):
        assert (mode == 'r'), "Can't convert TextIOBase to mode='rb'."
        yield file
        return

    if isinstance(file, io.IOBase):
        assert isinstance(file, (io.BufferedIOBase, io.TextIOBase))

        with contextlib.ExitStack() as stack:
            import gzip
            try:
                with seek_then_rewind(file, seek=0):
                    with gzip.open(file) as test:
                        test.read(2)
            except gzip.BadGzipFile:
                # Assume this is not a gzipped stream
                pass
            else:
                # gzip didn't complain
                # interpret this as a gzipped stream
                file = stack.enter_context(gzip.open(file))

            if (mode == 'r'):
                if not isinstance(file, io.TextIOBase):
                    assert isinstance(file, io.BufferedIOBase)
                    file = stack.enter_context(io.TextIOWrapper(file))

            yield file

        return

    if isinstance(file, str) and re.match(r"^(http|https|ftp)://", file):
        from idiva.download import download
        with download(file).now.open(mode='rb') as fd:
            with open_maybe_gz(fd, mode=mode) as fd:
                yield fd
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
