# RA, 2020-10-08

import os
import sys
from collections import defaultdict
from aligner03.utils import relpath, unlist1

import numpy as np
import typing

from unittest import TestCase
from pathlib import Path
from inclusive import range

from warnings import warn as warning



class TestPipeline(TestCase):
    def test_index(self):
        from aligner03.index import FM_Index as GenomeIndex

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
        GenomeIndex("N")
        GenomeIndex("ACTG")
        GenomeIndex("ACTGN")

        # Check the type
        loc_in_genome = GenomeIndex("ATTTTATTTG").query("TTT")
        self.assertIsInstance(loc_in_genome, list)
        for __ in loc_in_genome:
            self.assertIsInstance(__, int)

        self.assertCountEqual(GenomeIndex("ATTTTATTTTG").query("TTT"), [1, 2, 6, 7])

    def test_sanity(self):
        # import builtins
        # def print(*args, **kwargs):
        #     return builtins.print(*args, file=sys.stderr, **kwargs)

        data_root = Path(__file__).parent / "data_for_tests"
        source_path = data_root / "data_small"

        rs: np.random.RandomState
        rs = np.random.RandomState(0)

        #
        # STEP 1: LOAD THE REFERENCE GENOME
        #

        genome_file_fasta = unlist1(source_path.glob("genome*.fa"))

        print(F"Loading {relpath(genome_file_fasta)}")

        try:
            from aligner03.io import Read
            from aligner03.io import from_fasta
        except:
            warning("Using biopython for FASTA read.")
            from Bio import SeqIO
            from_fasta = (lambda f: SeqIO.read(f, format='fasta'))

        genome = from_fasta(genome_file_fasta)
        genome_seq = str(genome.seq)
        del genome

        EXPECTED_LENGTH = 4980
        self.assertIsInstance(genome_seq, str)
        self.assertEqual(len(genome_seq), EXPECTED_LENGTH)
        # all_reads: List[Read]

        #
        # STEP 2: INDEX THE REFERENCE GENOME
        #

        from aligner03.index import FM_Index as GenomeIndex

        try:
            GenomeIndex("ATTTTATTTG").query("T")
        except:
            warning(F"Basic genome index query fails")

        # Expect 0- based indexing
        self.assertListEqual(GenomeIndex("ATTTTT").query("ATT"), [0])

        GenomeIndex("ATTTTATTTTG").query("TTT")

        genome_index = GenomeIndex(genome_seq)

        query_the_genome = genome_index.query

        expected_read_length = 125
        expected_query_length = 20

        # HACK
        try:
            assert not query_the_genome(
                genome_seq[0:expected_query_length] +
                "".join(rs.choice(list("ACTG"), expected_query_length))
            )
        except AssertionError:
            __genome_index = GenomeIndex(genome_seq, k=expected_query_length)

            def query_the_genome(query):
                assert len(query) == expected_query_length
                return __genome_index.query(query, expected_query_length)
        finally:
            del genome_index

        self.assertEqual(query_the_genome(genome_seq[0:expected_query_length]), [0])

        #
        # STEP 3: LOAD SOME READS
        #

        read_file_fastq_12 = sorted(source_path.glob("*.fq"))
        self.assertEqual(len(read_file_fastq_12), 2)

        from aligner03.io import Read
        from_fastq: typing.Callable[[typing.AnyStr], typing.Iterable[Read]]

        try:
            from aligner03.io import from_fastq
            for read in from_fastq(read_file_fastq_12[0]):
                read.name
                read.desc
                read.seq
                read.phred
                self.assertEqual(len(read.seq), len(read.phred))
        except:
            raise
            warning("Using biopython for FASTQ file")
            from Bio import SeqIO
            def from_fastq(file):
                for read in SeqIO.parse(str(file), format='fastq'):
                    yield Read(read.name, read.description, read.seq, read.letter_annotations['phred_quality'])

        reads_by_file = {
            file.name: list(from_fastq(file))
            for file in read_file_fastq_12
        }

        for (f, reads) in reads_by_file.items():
            self.assertNotEqual(len(reads), 0)

        example_read: Read
        example_read = min(reads_by_file[min(reads_by_file.keys())], key=(lambda read: read.seq))


        try:
            # It's unlikely that the whole read is in the genome
            from aligner03.utils.strings import reverse, forward
            self.assertListEqual(query_the_genome(forward(example_read.seq)), [])
            self.assertListEqual(query_the_genome(reverse(example_read.seq)), [])
        except:
            warning("Either the whole read is there or another error.")

        # kmer length
        k = expected_query_length

        # https://en.wikipedia.org/wiki/Phred_quality_score
        #
        # | Phred Quality Score | Probability of incorrect base call | Base call accuracy |
        # |---------------------|------------------------------------|--------------------|
        # | 10                  | 1 in 10                            | 90%                |
        # | 20                  | 1 in 100                           | 99%                |
        # | 30                  | 1 in 1000                          | 99.9%              |
        # | 40                  | 1 in 10,000                        | 99.99%             |
        # | 50                  | 1 in 100,000                       | 99.999%            |
        # | 60                  | 1 in 1,000,000                     | 99.9999%           |

        def all_kmers(read: Read) -> typing.Dict[float, typing.List[tuple]]:
            """
            Get all kmers arranged by phred score.
            """
            by_score = defaultdict(list)
            for i in range[0, len(read.seq) - k]:
                ii = slice(i, i + k)
                by_score[
                    np.average(read.phred[ii])
                ].append(
                    (i, read.seq[ii])
                )
            return dict(by_score)

        # print("These are all kmers from the read:", *sorted(all_kmers(example_read).items()), sep="\n > ")

        def propose_mapping(kmers_by_phred: dict):
            for (phred, kmers) in sorted(kmers_by_phred.items(), reverse=True):
                for (i, kmer) in kmers:
                    for j in query_the_genome(kmer):
                        yield (phred, kmer, i, j)

        proposed_mappings = [
            (read, list(propose_mapping(all_kmers(read))))
            for read in [example_read, example_read.reverse]
        ]

        for (read, mappings) in proposed_mappings:
            print(F"Number of proposed mappings (forward = {read.is_forward}):", len(mappings))

        # Choose which to go with
        # Note: overwrites these variables
        (example_read, proposed_mappings) = max(
            proposed_mappings,
            key=(lambda xx: len(xx[1]))
        )


        for (i, kmer, phred, j) in proposed_mappings:
            # print(F"{kmer} at {i} and quality {phred} maps to reference:", j)
            self.assertEqual(genome_seq[j:(j + len(kmer))], kmer)

        # Take any mapping
        (phred, kmer, i, j) = min(proposed_mappings)

        # Perform local alignment of the read in the vicinity of j
        grace_margin = expected_read_length // 2
        jj = slice(j - i - grace_margin, j - i + expected_read_length + grace_margin)
        # Note: this fails if j is close to the boundary of the genome
        self.assertEqual(len(genome_seq[jj]), expected_read_length + grace_margin * 2)


        #
        # STEP 4: GET AN ALIGNMENT
        #

        from aligner03.align import SmithWaterman as Aligner
        aligner = Aligner()

        print(F"Trying to align {example_read.seq} within {genome_seq[jj]}")

        alignments = list(aligner(genome_seq[jj], example_read.seq))

        for alignment in alignments:
            print(" >", alignment)


        #
        # STEP 5: DUMP TO SAM
        #

        # Unfinished


        # TODO
        # Semi-local alignment rather than local Smith-Waterman
        # The cognate reads should be nearby
        #
