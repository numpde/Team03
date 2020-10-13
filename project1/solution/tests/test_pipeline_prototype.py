# RA, 2020-10-08

import os
import sys
from collections import defaultdict
from humdum.utils import relpath, unlist1

import numpy as np
import typing

from unittest import TestCase
from pathlib import Path
from inclusive import range

from warnings import warn as warning



class TestPipeline(TestCase):
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

        from humdum.io import from_fasta
        genome_seq = unlist1(from_fasta(genome_file_fasta)).seq

        EXPECTED_LENGTH = 4980
        self.assertIsInstance(genome_seq, str)
        self.assertEqual(len(genome_seq), EXPECTED_LENGTH)

        #
        # STEP 2: INDEX THE REFERENCE GENOME
        #

        from humdum.index import FmIndex as GenomeIndex

        try:
            GenomeIndex("ATTTTATTTG").query("T")
        except:
            warning(F"Basic genome index query fails")

        # Expect 0- based indexing
        self.assertListEqual(GenomeIndex("ATTTTT").query("ATT"), [0])

        GenomeIndex("ATTTTATTTTG").query("TTT")

        query_the_genome = GenomeIndex(genome_seq).query

        expected_read_length = 125
        query_length = 26

        self.assertEqual(query_the_genome(genome_seq[0:query_length]), [0])

        #
        # STEP 3: LOAD SOME READS
        #

        read_file_fastq_12 = sorted(source_path.glob("*.fq"))
        self.assertEqual(len(read_file_fastq_12), 2)

        from humdum.io import Read
        from_fastq: typing.Callable[[typing.AnyStr], typing.Iterable[Read]]

        from humdum.io import from_fastq
        for read in from_fastq(read_file_fastq_12[0]):
            # print(read.name, read.desc, read.seq, read.phred)
            self.assertEqual(len(read.seq), len(read.phred))

        reads_by_file = {
            file.name: list(from_fastq(file))
            for file in read_file_fastq_12
        }

        # Check that we got some reads from each file
        for (f, reads) in reads_by_file.items():
            self.assertNotEqual(len(reads), 0)

        example_read: Read
        example_read = min(reads_by_file[min(reads_by_file.keys())], key=(lambda read: read.seq))


        try:
            # It's unlikely that the whole read is in the genome
            from humdum.utils.strings import reverse, forward
            self.assertListEqual(query_the_genome(forward(example_read.seq)), [])
            self.assertListEqual(query_the_genome(reverse(example_read.seq)), [])
        except:
            warning("Either the whole read is there or another error.")

        # kmer length
        k = query_length

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

        from humdum.map import all_kmers_by_score

        # print(
        #     "These are all kmers from the read:",
        #     *sorted(all_kmers_by_score(example_read).items()),
        #     sep="\n > "
        # )

        def propose_mapping(kmers_by_phred: dict):
            for (phred, kmers) in sorted(kmers_by_phred.items(), reverse=True):
                for (i, kmer) in kmers:
                    for j in query_the_genome(kmer):
                        yield (phred, kmer, i, j)

        proposed_mappings = [
            (read, list(propose_mapping(all_kmers_by_score(read, k))))
            for read in [example_read, example_read.reversed]
        ]

        for (read, mappings) in proposed_mappings:
            print(F"Number of proposed mappings ({'Forward' if read.is_forward else 'Backward'}):", len(mappings))

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
        from humdum.map import propose_window
        (a, b) = propose_window(
            read_length=len(example_read.seq),
            read_loc=i,
            ref_length=len(genome_seq),
            ref_loc=j,
        )

        # Slices into the reference genome
        jj = slice(a, b)

        print("Range on the genome to align the read: ", [a, b])

        #
        # STEP 4: GET AN ALIGNMENT
        #

        from humdum.align import SmithWaterman as Aligner
        aligner = Aligner()

        print(F"Trying to align {example_read.seq} within {genome_seq[jj]}")

        alignments = list(aligner(ref=genome_seq[jj], query=example_read.seq))

        for alignment in alignments:
            print(" >", alignment)


        #
        # STEP 5: CHECK WITH SAM
        #

        from pysam import AlignedSegment
        from humdum.io import from_sam_pysam
        from humdum.io.sam import Flag
        segment: AlignedSegment
        for segment in from_sam_pysam(unlist1(source_path.glob("*.sam"))):
            # print(segment.query_name, segment.flag, segment.qname)
            flag = Flag(segment.flag)
            name = segment.query_name + {True: "/1", False: "/2"}[flag.is_minus_strand]
            if (example_read.name == name):
                print("Reference alignment:", segment, sep='\n')


        # TODO
        # Semi-local alignment rather than local Smith-Waterman
        # The cognate reads should be nearby
        #
