# JG, 2020-10

from aligner03.io.fasta__orphaned import Fasta
import filecmp
import io
import sys
import unittest

class TestFasta(unittest.TestCase):
#all fa tests
    def test_missingFile(self):
    #test for not existing file
        with self.assertRaises(SystemExit) as cm:
            Fasta("test123.txt")
        self.assertEqual(cm.exception.code, "Could not open the following file: test123.txt")
    def test_non_positive(self):
    #test non-positive seq-lengths
        for i in ["hi",0,-1]:
            fasta=Fasta("tests\data\empty.txt")
            with self.assertRaises(SystemExit) as cm:
                for j in fasta.yield_sequences(i):
                    pass
            self.assertEqual(cm.exception.code,"Sequence length must be positive integer.")
    def test_empty(self):
    # test on empty file
        fasta=Fasta("tests\data\empty.txt")
        self.assertEqual(fasta.title,"")
        variable="some testcase"
        for i in fasta.yield_sequences(1):
            variable="another testcase"
        self.assertEqual(variable,"some testcase")
    def test_too_small(self):
    # test if sequence is too small to get a long enough sequence
        fasta=Fasta('tests\data\io\\fasta_too_small.txt')
        self.assertEqual(fasta.title,"asdf")
        variable="some testcase"
        for i in fasta.yield_sequences(4):
            variable="another testcase"
        self.assertEqual(variable,"some testcase")
    def test_small_data(self):
    # test on small data
        fasta=Fasta("tests\data\io\genome.chr22.5K.fa")
        self.assertEqual(fasta.title,">22_5K")
        with open("tests\data\io\genome_correct.txt","r",encoding="utf-8") as correct:
            for i in fasta.yield_sequences(4):
                self.assertEqual(str(i),correct.readline().rstrip())
if __name__ == '__main__':
    unittest.main()