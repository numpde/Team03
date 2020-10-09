# RA, 2020-10-09

from aligner03.index.fm import FM_Index

from unittest import TestCase


class TestIndex(TestCase):
    def test_constructor(self):
        pass

    def test_dollar(self):
        pass

    def test_perfect_match_mini(self):
        ref = "TAGAGAGATCGATTTTTTCTTGACTGACTGACTCAG"

        for query in ["ACT", "T", "TTT", "TTTT"]:
            fm_index = FM_Index(ref)
            hits = fm_index.query(query)

            # Test for precision
            for i in hits:
                self.assertEqual(ref[i:(i + len(query))], query)

            # Test for recall
            from aligner03.utils import find_all
            self.assertCountEqual(hits, find_all(template=ref, pattern=query))
