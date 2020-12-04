# RA, 2020-10-13

import io

from unittest import TestCase
from pathlib import Path
from tempfile import NamedTemporaryFile

from idiva.io import open_maybe_gz

data_root = Path(__file__).parent / "data_for_tests/gz"

ref_file = data_root / "test.txt"
ref_file_gz = data_root / "test.txt.gz"

with open(ref_file, mode='r') as fd_ref:
    reference = fd_ref.read()


class TestGz(TestCase):
    def test_plain(self):
        with open_maybe_gz(ref_file, mode='r') as fd:
            self.assertEqual(reference, fd.read())

    def test_gz(self):
        with open_maybe_gz(ref_file_gz, mode='r') as fd:
            candidate = fd.read()
            self.assertIsInstance(candidate, str)

        self.assertEqual(reference, candidate)

    def test_open_rb(self):
        with open_maybe_gz(ref_file, mode='rb') as fd:
            with io.TextIOWrapper(fd) as fd:
                self.assertEqual(reference, fd.read())

    def test_dont_open_r_as_rb(self):
        with open(ref_file, mode='r') as fd:
            with self.assertRaises(AssertionError):
                with open_maybe_gz(fd, mode='rb'):
                    pass

    def test_open_from_rb_as_r(self):
        with open(ref_file, mode='rb') as fd:
            with open_maybe_gz(fd, mode='r') as fd:
                self.assertEqual(reference, fd.read())

    def test_open_from_rb_as_rb(self):
        with open(ref_file, mode='rb') as fd:
            with open_maybe_gz(fd, mode='rb') as fd:
                with io.TextIOWrapper(fd) as fd:
                    self.assertEqual(reference, fd.read())
