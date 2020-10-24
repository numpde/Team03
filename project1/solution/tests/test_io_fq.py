# RA, 2020-10-09

import unittest
from pathlib import Path

from humdum.io import from_fastq, Read
from humdum.utils import relpath

data_root = Path(__file__).parent / "data_for_tests"


class TestFastqReader(unittest.TestCase):
    def test_data_small_vs_biopython(self):
        source_path = data_root / "data_small"

        from Bio import SeqIO

        files = source_path.glob("*.fq")

        for file in files:
            reference_reads = list(SeqIO.parse(file, format='fastq'))
            candidate_reads = list(from_fastq(file))
            self.assertEqual(len(reference_reads), len(candidate_reads))
            for (reference, candidate) in zip(reference_reads, candidate_reads):
                self.assertIsInstance(reference, SeqIO.SeqRecord)
                self.assertIsInstance(candidate, Read)
                self.assertEqual(str(reference.seq), candidate.seq)
                self.assertEqual(reference.name, candidate.name)
                self.assertEqual(reference.letter_annotations["phred_quality"], candidate.phred)
                # self.assertEqual(reference.description, candidate.desc)

    def test_on_data_large(self):
        source_path = data_root / "data"

        files = source_path.glob("*.fq*")
        for file in files:
            n = len(list(from_fastq(file)))
            print(F"File {relpath(file)} has {n} reads.")
