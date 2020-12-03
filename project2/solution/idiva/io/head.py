# RA, 2020-12-02

import io
import contextlib


@contextlib.contextmanager
def head(src: io.TextIOBase, n=100) -> io.TextIOBase:
    from idiva.utils import seek_then_rewind
    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(mode='w') as tf:
        with seek_then_rewind(src):
            # Write meta and header

            line = "##"
            while line.startswith("##"):
                line = src.readline().strip()
                print(line, file=tf)

            assert line.startswith("#")

            # Write a few datalines

            for _ in range(n):
                line = src.readline().strip()
                assert not line.startswith("#")
                print(line, file=tf)

        with open(tf.name, mode='r') as fd:
            yield fd
