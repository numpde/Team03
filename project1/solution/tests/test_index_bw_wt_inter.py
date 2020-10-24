
from unittest import TestCase

# Note: there could be an issue with multiple TestCase classes per file [RA]
from humdum.index.bw import BurrowsWheeler
from humdum.index.wt import WaveletTree


class TestIndexBWWT(TestCase):
    def test_constructor(self):
        pass

    def test_rank(self):
        ref_genome = "ACGCATGATCACTAGCTAGCATCGACCNANCN"

        bw = BurrowsWheeler(ref_genome)
        wt = WaveletTree(ref_genome)

        for char in ['A', 'C', 'G', 'T', 'N']:
            for i in range(len(ref_genome)):
                self.assertEqual(bw.rank(char, i), wt.rank(char, i))

    def test_sa(self):
        ref_genome = "GTCAACGCATGATCGATACGCATGATCGACCNANCN"

        bw = BurrowsWheeler(ref_genome)
        wt = WaveletTree(ref_genome)

        for i in range(len(ref_genome)):
            self.assertEqual(bw.get_sa(i), wt.get_sa(i))
