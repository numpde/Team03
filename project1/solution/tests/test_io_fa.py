# RA, 2020-10-09

import unittest
from pathlib import Path

from humdum.io import from_fasta, Sequence

data_root = Path(__file__).parent / "data_for_tests"


class TestFastaReader(unittest.TestCase):
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
        files = source_path.glob("*.fa*")
        assert files

        for file in files:
            for genome in from_fasta(file):
                print(F"Got genome of length {len(genome.seq)}")
