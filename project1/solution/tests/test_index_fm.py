# RA, 2020-10-09

from humdum.index import FmIndex as GenomeIndex

from unittest import TestCase


class TestIndex(TestCase):
    def test_constructor(self):
        GenomeIndex("ACGT")

    def test_returns_original(self):
        original = "ACGT"
        self.assertEqual(original, str(GenomeIndex(original)))

        original = "ACGT" * 100
        self.assertEqual(original, str(GenomeIndex(original)))

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

        self.assertListEqual(GenomeIndex("ATTTTATTTG").query("TTTTTT"), [])

        GenomeIndex("A")
        GenomeIndex("ACTG")

        # Check the type
        loc_in_genome = GenomeIndex("ATTTTATTTG").query("TTT")
        self.assertIsInstance(loc_in_genome, list)
        for __ in loc_in_genome:
            self.assertIsInstance(__, int)

        self.assertCountEqual(GenomeIndex("ATTTTATTTTG").query("TTT"), [1, 2, 6, 7])

        with self.assertRaises(NotImplementedError):
            GenomeIndex("N")

        with self.assertRaises(NotImplementedError):
            GenomeIndex("ACTGN")

    def test_perfect_match_mini(self):
        ref = "TAGAGAGATCGATTTTTTCTTGACTGACTGACTCAG"

        for query in ["ACT", "T", "TTT", "TTTT"]:
            fm_index = GenomeIndex(ref)
            hits = fm_index.query(query)

            # Test for precision
            for i in hits:
                self.assertEqual(ref[i:(i + len(query))], query)

            # Test for recall
            from humdum.utils import find_all
            self.assertCountEqual(hits, find_all(template=ref, pattern=query))

    def test_read_write(self):

        ref = "TAGAGAGACGTACGATCGACTGACTAGCTAGCACACACACTGACTCGATCGAC"

        fm_index = GenomeIndex(ref)

        fm_index.write("data_for_tests/index_data/")
        fm_index2 = GenomeIndex.read("data_for_tests/index_data/")

        self.assertIsInstance(fm_index2,GenomeIndex)

        self.assertEqual(fm_index.bwt.code, fm_index2.bwt.code)
        self.assertEqual(fm_index.bwt.sa, fm_index2.bwt.sa)
        self.assertEqual(fm_index.bwt.decode(), fm_index2.bwt.decode())

        self.assertEqual(fm_index.F, fm_index2.F)
        self.assertEqual(fm_index.tally, fm_index2.tally)
        self.assertEqual(fm_index.next_chars._data, fm_index2.next_chars._data)

