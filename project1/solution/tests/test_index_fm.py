
from humdum.index import FmIndex as GenomeIndex

from unittest import TestCase
from pathlib import Path

data_root = Path(__file__).parent / "data_for_tests"


class TestIndex(TestCase):
    def test_constructor(self):
        GenomeIndex("ACGT")

    def test_basic(self):

        with self.assertRaises(ValueError):
            GenomeIndex("ACTG$")

        with self.assertRaises(ValueError):
            GenomeIndex("$ACTG")

        with self.assertRaises(ValueError):
            GenomeIndex("#ACTG")

        with self.assertRaises(ValueError):
            GenomeIndex("")

        with self.assertRaises(ValueError):
            GenomeIndex("ACGT").query("ACGTTTTTTTTT")

        with self.assertRaises(ValueError):
            GenomeIndex("ACGT").query("GT$")

        with self.assertRaises(ValueError):
            GenomeIndex("A", -1)

        with self.assertRaises(ValueError):
            GenomeIndex("A", compression_occ=-1)

        with self.assertRaises(ValueError):
            GenomeIndex("A", compression_sa=-1)

        self.assertListEqual(GenomeIndex("ATTTTATTTG").query("TTTTTT"), [])

        GenomeIndex("A")
        GenomeIndex("ACTG")

        # Check the type
        loc_in_genome = GenomeIndex("NANNTTTTATTTG").query("TTT")
        self.assertIsInstance(loc_in_genome, list)
        for __ in loc_in_genome:
            self.assertIsInstance(__, int)

        self.assertCountEqual(GenomeIndex("ATTTTANNTTTTGN").query("TTT"), [1, 2, 8, 9])

    def test_perfect_match_mini(self):
        ref = "TAGANNGAGATCNGNNNATTNTTTNTCTNNNTGANCTGNACTGACTCAAAAAGN"

        for query in ["ACT", "T", "TTT", "TTTT", "GNACT", "AAAA"]:
            fm_index = GenomeIndex(ref)
            hits = fm_index.query(query)

            # Test for precision
            for i in hits:
                self.assertEqual(ref[i:(i + len(query))], query)

            # Test for recall
            from humdum.utils import find_all
            self.assertCountEqual(hits, find_all(template=ref, pattern=query))

    def test_match_with_compression(self):
        ref = "TAGAATCGTTTTTTTTTTATCGACTACNACTACAAAAAAAAATGATCNTACNGTAANNNNNTTNTTTNTCTNNNTGANCTGNACTGACTCAAAAAGN"

        for comp in range(1, 100):

            fm_index = GenomeIndex(ref, compression_occ=comp, compression_sa=comp + 1)
            for query in ["ACT", "T", "TTT", "TTTT", "GNACT", "AAAA"]:
                hits = fm_index.query(query)

                # Test for precision
                for i in hits:
                    self.assertEqual(ref[i:(i + len(query))], query)

                # Test for recall
                from humdum.utils import find_all
                self.assertCountEqual(hits, find_all(template=ref, pattern=query))

    def test_decoding(self):
        refs = ["TAGAATCGTTTTTTTTTTATCGACTACNACTACAAAAAAAAATGATCNTACNGTAA",
                "TTTTTTTTTTTAAAAAAAAAACCCCCCCGGN",
                "AGCTA", "T"]
        for ref in refs:
            for comp in range(1, 50):
                fm_index = GenomeIndex(ref, compression_occ=comp, compression_sa=comp)
                self.assertEqual(ref, str(fm_index))

    def test_read_write(self):

        ref = "NTAGAGNANGACGTACNGATCGANCTGACTNAGCTNNNAGCACACACACTGACTCNNGATCGACNNN"

        fm_index = GenomeIndex(ref)

        Path(data_root / "index_data").mkdir(parents=True, exist_ok=True)

        fm_index.write(data_root / "index_data/index_small.data")
        fm_index2 = GenomeIndex.read(data_root / "index_data/index_small.data")

        self.assertIsInstance(fm_index2, GenomeIndex)

        self.assertEqual(fm_index.bwt.sa, fm_index2.bwt.sa)
        self.assertEqual(str(fm_index), str(fm_index2))
        self.assertEqual(fm_index.bwt.f, fm_index2.bwt.f)
        self.assertEqual(fm_index.bwt.next_chars._data, fm_index2.bwt.next_chars._data)

    def test_wavelet(self):

        ref = "TAGANNGAGATCNGNNNATTNTTTNTCTNNNTGANCTGNACTGACTCAAAAAGN"

        for query in ["ACT", "T", "TTT", "TTTT", "GNACT", "AAAA"]:
            fm_index = GenomeIndex(ref)
            wavelet = GenomeIndex(ref, wavelet=True)
            fm_hits = fm_index.query(query)
            wavelet_hits = wavelet.query(query)
            self.assertEqual(fm_hits, wavelet_hits)
