
import unittest


class TestIndex(unittest.TestCase):
    def test_imports_naive(self):
        from humdum.index import NaiveIndex as GenomeIndex
        self.assertIsNotNone(GenomeIndex)

    def test_imports_fancy(self):
        from humdum.index import FmIndex as GenomeIndex
        self.assertIsNotNone(GenomeIndex)

    def test_healthcheck_naive(self):
        from humdum.index import NaiveIndex as GenomeIndex
        from humdum.qc import healthcheck_index
        healthcheck_index(GenomeIndex)

    def test_healthcheck_fancy(self):
        from humdum.index import FmIndex as GenomeIndex
        from humdum.qc import healthcheck_index
        healthcheck_index(GenomeIndex)
