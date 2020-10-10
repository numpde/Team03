# RA, 2020-10-09

import unittest
from pathlib import Path

from humdum.io import from_fasta, Sequence


class TestFastaReader(unittest.TestCase):
    def test_on_data_small(self):
        data_root = Path(__file__).parent / "data_for_tests"
        source_path = data_root / "data_small"

        from Bio import SeqIO

        files = source_path.glob("*.fa")

        for file in files:
            reference_reads = list(SeqIO.parse(file, format='fasta'))
            candidate_reads = list(from_fasta(file))
            self.assertEqual(len(reference_reads), len(candidate_reads))
            for (reference, candidate) in zip(reference_reads, candidate_reads):
                self.assertIsInstance(reference, SeqIO.SeqRecord)
                self.assertIsInstance(candidate, Sequence)
                self.assertEqual(str(reference.seq), candidate.seq)
                self.assertEqual(reference.description, candidate.desc)
