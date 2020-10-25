# RA, 2020-10-25

import unittest


class TestIndex(unittest.TestCase):
    def test_imports_naive(self):
        from humdum.index import NaiveIndex as GenomeIndex
        self.assertIsNotNone(GenomeIndex)

    def test_imports_fancy(self):
        from humdum.index import FmIndex as GenomeIndex
        self.assertIsNotNone(GenomeIndex)
