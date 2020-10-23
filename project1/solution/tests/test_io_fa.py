# RA, 2020-10-09

import unittest
from pathlib import Path
from tempfile import NamedTemporaryFile

from humdum.io import from_fasta, Sequence
from humdum.utils import first

data_root = Path(__file__).parent / "data_for_tests"


class TestFastaReader(unittest.TestCase):
    def test_reads_well(self):
        desc = "Hello"
        seq1 = "ABC"
        seq2 = "DEF"
        with NamedTemporaryFile(mode='w') as fn:
            print(*[">" + desc, seq1, seq2], sep='\n', file=fn, flush=True)
            record = first(from_fasta(fn.name))
            self.assertEqual(record.desc, desc)
            self.assertEqual(record.seq, seq1 + seq2)

    def test_fails_when_many(self):
        with NamedTemporaryFile(mode='w') as fn:
            print(*[">A", "N", ">B", "N"], sep='\n', file=fn, flush=True)
            with self.assertRaises(AssertionError):
                list(from_fasta(fn.name))

    def test_on_data_small(self):
        source_path = data_root / "data_small"
        files = list(source_path.glob("*.fa"))
        assert files

        from Bio import SeqIO

        for file in files:
            reference_reads = list(SeqIO.parse(file, format='fasta'))
            candidate_reads = list(from_fasta(file))
            self.assertEqual(len(reference_reads), len(candidate_reads))
            for (reference, candidate) in zip(reference_reads, candidate_reads):
                self.assertIsInstance(reference, SeqIO.SeqRecord)
                self.assertIsInstance(candidate, Sequence)
                self.assertEqual(str(reference.seq), candidate.seq)
                self.assertEqual(reference.description, candidate.desc)

    def test_on_data_big(self):
        source_path = data_root / "data"
        files = list(source_path.glob("*.fa.gz"))
        assert files

        for file in files:
            for genome in from_fasta(file):
                self.assertEqual(len(genome.seq), 51304566)
                self.assertTrue(genome.seq.strip("N").endswith("CGGATT"))
