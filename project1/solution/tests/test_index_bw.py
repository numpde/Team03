
from unittest import TestCase

# Note: there could be an issue with multiple TestCase classes per file [RA]
from humdum.index.bw import BurrowsWheeler


class TestIndexBW(TestCase):
    def test_constructor(self):
        pass

    def test_basic(self):
        with self.assertRaises(ValueError):
            BurrowsWheeler("")

        with self.assertRaises(ValueError):
            BurrowsWheeler("A", "A")

        with self.assertRaises(ValueError):
            BurrowsWheeler("A", strategy="fun")

        with self.assertRaises(ValueError):
            BurrowsWheeler("A", compression_occ=-1)

        with self.assertRaises(ValueError):
            BurrowsWheeler("A", compression_sa=-1)

        BurrowsWheeler("CACGTACGTGTGCTAACACGTGTGTTTTTGAC")

        suffix = BurrowsWheeler("GCAGTN").sa
        encoded = BurrowsWheeler("ACGTGTAC").get_bwt("ACGTGTAC")

        self.assertIsInstance(encoded, str)
        self.assertIsInstance(suffix, list)

        string = "ACGATCGATCAGTAC"
        self.assertEqual(len(string), len(BurrowsWheeler(string)))

    def test_encode_decode(self):
        string = "NNCACGTACGTGTGCTAACACGTGTGTTTTTGAC"
        bwt = BurrowsWheeler(string)
        self.assertEqual(str(bwt), string)

    def test_sa_compression(self):

        refs = ["TAGAATCGTTTTTTTTTTATCGACTACNACTACAAAAAAAAATGATCNTACNGTAA",
                "TTTTTTTTTTTAAAAAAAAAACCCCCCCGGN",
                "AGCTA", "T"]
        for ref in refs:
            for comp in range(1, 50):
                for i in range(1, len(ref)):
                    bw_uncompressed = BurrowsWheeler(ref, compression_occ=1, compression_sa=1)
                    bw = BurrowsWheeler(ref, compression_occ=1, compression_sa=comp)
                    self.assertEqual(bw_uncompressed.sa[i], bw.get_sa(i))

    def test_algorithms(self):

        ref = "TAGAATCGTTTTTTTTTTATCGACTACNACTACAAAAAAAAATGATCNTACNGTAATTTTTTTTTTTAAAAAAAAAACCCCCCCGGN"

        simple = BurrowsWheeler(ref, strategy="Simple")
        manber = BurrowsWheeler(ref, strategy="ManberMyers")
        kaerkkaeinen = BurrowsWheeler(ref, strategy="KaerkkaeinenSanders")

        self.assertEqual(str(simple), str(manber))
        self.assertEqual(str(simple), str(kaerkkaeinen))
