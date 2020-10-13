# RA, 2020-10-13

from unittest import TestCase
from pathlib import Path

from humdum.io import open_maybe_gz

data_root = Path(__file__).parent / "data_for_tests/gz"


class TestGz(TestCase):
    def test_gz(self):
        with open(data_root / "test.txt", mode='r') as fd:
            reference = fd.read()

        with open_maybe_gz(data_root / "test.txt.gz") as fd:
            candidate = fd.read()
            self.assertIsInstance(candidate, str)

        self.assertEqual(reference, candidate)

