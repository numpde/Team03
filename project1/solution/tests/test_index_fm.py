# RA, 2020-10-09

from aligner03.index import FmIndex as GenomeIndex

from unittest import TestCase


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
            fm_index = FmIndex(ref)
            hits = fm_index.query(query)

            # Test for precision
            for i in hits:
                self.assertEqual(ref[i:(i + len(query))], query)

            # Test for recall
            from aligner03.utils import find_all
            self.assertCountEqual(hits, find_all(template=ref, pattern=query))
