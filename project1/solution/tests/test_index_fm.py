
from aligner03.index.fm import FM_Index

from unittest import TestCase


class TestIndex(TestCase):
    def test_constructor(self):
        pass

    def test_dollar(self):
        pass

    def test_perfect_match(self):
        ref = "TAGAGAGATCGATCGACTGACTGACTCAG"
        query = "ACT"

        fm_index = FM_Index(ref)
        hits = fm_index.query(query)

        for i in hits:
            self.assertEqual(ref[i:(i + len(query))], query)

        def find(pattern, template):
            """
            Find all occurrences of pattern string in template string.
            Note: returns a list of starting positions with 1- based indexing.
            """
            import re
            return [(m.start() + 1) for m in re.finditer("(?=" + pattern + ")", template)]


        # TODO:
        # Make sure that _all_ occurrences are returned

