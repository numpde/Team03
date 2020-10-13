# RA, 2020-10-13

from unittest import TestCase
from humdum.io.sam import Flag


class TestFlag(TestCase):
    def test_switch_is_minus_strand(self):
        f = Flag(0)
        for p in [0, 1, 1, 0, 1, 1, 0, 0, 0, 1]:
            f.is_minus_strand = bool(p)
            self.assertEqual(p, f.is_minus_strand)
