
from unittest import TestCase

# Note: there could be an issue with multiple TestCase classes per file [RA]
from humdum.index.bw import BurrowsWheeler
from typing import List


class TestIndexBW(TestCase):
    def test_constructor(self):
        pass

    def test_basic(self):
        with self.assertRaises(ValueError):
            BurrowsWheeler("")

        with self.assertRaises(ValueError):
            BurrowsWheeler("A", "A")

        BurrowsWheeler("some string")
        BurrowsWheeler("some string with special characters .:,;{$£}[¨!]+*ç%&/()=?|@#¼½¬|")

        encoded = BurrowsWheeler("some string").code
        suffix = BurrowsWheeler("some string").sa

        self.assertIsInstance(encoded, str)
        self.assertIsInstance(suffix, list)

    def test_encode_decode(self):

        self.assertEqual(BurrowsWheeler("simple string").decode(),"simple string")
        self.assertEqual(BurrowsWheeler("some string with special characters!+*/-.:,;-_+*ç%&/()=?è!£àé<>\}][|@#¼½¬|¢~").decode(),
                         "some string with special characters!+*/-.:,;-_+*ç%&/()=?è!£àé<>\}][|@#¼½¬|¢~")
        self.assertEqual(BurrowsWheeler("This is a longer text to test Burrows Wheeler Transform. "
                                        "Wikipedia: The Burrows–Wheeler transform (BWT, also called block-sorting compression) "
                                        "rearranges a character string into runs of similar characters. "
                                        "This is useful for compression, since it tends to be easy to compress a string"
                                        " that has runs of repeated characters by techniques such as move-to-front"
                                        " transform and run-length encoding. More importantly, the transformation"
                                        " is reversible, without needing to store any additional data except"
                                        " the position of the first original character. The BWT is thus a free "
                                        "method of improving the efficiency of text compression algorithms,"
                                        " costing only some extra computation. ").decode(),
                         "This is a longer text to test Burrows Wheeler Transform. "
                         "Wikipedia: The Burrows–Wheeler transform (BWT, also called block-sorting compression) "
                         "rearranges a character string into runs of similar characters. "
                         "This is useful for compression, since it tends to be easy to compress a string"
                         " that has runs of repeated characters by techniques such as move-to-front"
                         " transform and run-length encoding. More importantly, the transformation"
                         " is reversible, without needing to store any additional data except"
                         " the position of the first original character. The BWT is thus a free "
                         "method of improving the efficiency of text compression algorithms,"
                         " costing only some extra computation. ")




