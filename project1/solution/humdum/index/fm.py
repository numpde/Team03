# LB, pre- 2020-10-09


# TODO:
# -BWT can be compressed => Wavelet trees (also occurrence matrix could be removed)

# DONE:
# -First column compressed to a dict of 6 entries
# -suffix compressed (additional bitvector)
# -occurrence matrix compressed (store every kth row)
import time
from sys import getsizeof

from humdum.index import BurrowsWheeler
from typing import List
import bz2
import pickle


"""
Paolo Ferragina and Giovanni Manzini
An experimental study of a compressed index_data
InformationSciences, 135(1-2):13–28, June 2001
ISSN 00200255
doi:10.1016/S0020-0255(01)00098-6
URL: http://linkinghub.elsevier.com/retrieve/pii/S0020025501000986.
"""

"""
How the FM-Index works by Ben Langmead
https://www.youtube.com/watch?v=kvVGj5V65io
https://www.youtube.com/watch?v=UHGgpfxlaiE
https://www.youtube.com/watch?v=DYyXgxYmSYU
"""


class FmIndex:
    """
    Implements the FM-index_data for substring lookup.

    Upon creation a Burrows Wheeler object is created

    The compression for the occurrence matrix (compression_occ)
    and the compression for the suffix array (compression_sa) can be select upon creation of the object

    The FM-index_data works only for strings over the alphabet {A, C, G, N, T}.
    """

    def __string_checks(self, string):
        if '$' in string:
            raise ValueError("The string should not contain '$'.")

        if not set(string).issubset(set("ACGTN")):
            raise ValueError("The string contains unexpected characters.")

        if not string:
            raise ValueError("The string may not be empty.")

    def __init__(self, reference_genome: str, compression_occ: int = 8, compression_sa: int = 1):

        if compression_occ < 1 or compression_sa < 1:
            raise ValueError("compression coefficients need to be strictly positive >=0")

        self.__string_checks(reference_genome)

        self.compression_occ = compression_occ

        reference_genome += "$"

        # contains the Suffix Array of the reference genome
        # can create the Burrows Wheeler Transformation on demand
        # size: len(ref_genome) * sizeof(int)
        self.bwt = BurrowsWheeler(reference_genome, compression_occ=compression_occ, compression_sa=compression_sa)

    def __len__(self):
        """
        The length of the reference string.
        """
        return len(self.bwt)

    def __str__(self):
        """
        returns the original string
        """
        return str(self.bwt)

    def query_new(self, sample: str) -> List[int]:
        """
                Query index_data for a given sample.
                Returns a list of all positions in the reference genome (0-based).
                """

        self.__string_checks(sample)

        if len(sample) >= len(self.bwt.code):
            raise ValueError(
                F"Sample length may not exceed that of the reference genome ({len(self)})."
            )

        half_compression = 0.5 * self.compression_occ

        n = len(self.bwt.code) - 1

        i = len(sample) - 1
        c = sample[i]
        sp = self.bwt.f[c] - 1
        ep = self.bwt.f[self.bwt.next_chars[c]] - 1

        div = 1. / self.compression_occ

        while ep > sp and i >= 1:
            print(i)
            c = sample[i - 1]

            if sp % self.compression_occ == 0:
                rank_sp = self.bwt.tally[c][int(sp * div)]

            elif sp % self.compression_occ < half_compression or sp > int(n * div) * self.compression_occ:
                count = 0
                for up in range(sp, sp - (sp % self.compression_occ), -1):
                    if self.bwt.code[up] == c:
                        count += 1

                rank_sp = self.bwt.tally[c][int(sp * div)] + count

            else:
                count = 0
                for down in range(sp + 1, sp + (self.compression_occ - (sp % self.compression_occ)) + 1):
                    if self.bwt.code[down] == c:
                        count += 1

                rank_sp = self.bwt.tally[c][int(sp * div + 1)] - count

            if ep % self.compression_occ == 0:
                rank_ep = self.bwt.tally[c][int(ep * div)]

            elif ep % self.compression_occ < half_compression or ep > int(n * div) * self.compression_occ:
                count = 0
                for up in range(ep, ep - (ep % self.compression_occ), -1):
                    if self.bwt.code[up] == c:
                        count += 1

                rank_ep = self.bwt.tally[c][int(ep * div)] + count

            else:
                count = 0
                for down in range(ep + 1, ep + (self.compression_occ - (ep % self.compression_occ)) + 1):
                    if self.bwt.code[down] == c:
                        count += 1
                rank_ep = self.bwt.tally[c][int(ep * div + 1)] - count

            sp = self.bwt.f[c] + rank_sp - 1
            ep = self.bwt.f[c] + rank_ep - 1
            i = i - 1

        if sp >= ep:
            return []
        else:
            return [self.bwt.get_sa(i) for i in range(sp + 1, ep + 1)]


    def query(self, sample: str) -> List[int]:
        """
        Query index_data for a given sample.
        Returns a list of all positions in the reference genome (0-based).
        """

        self.__string_checks(sample)

        if len(sample) >= len(self.bwt.code):
            raise ValueError(
                F"Sample length may not exceed that of the reference genome ({len(self)})."
            )

        half_compression = 0.5 * self.compression_occ
        n =  len(self.bwt.code)- 1

        i = len(sample) - 1
        c = sample[i]
        sp = self.bwt.f[c] - 1
        ep = self.bwt.f[self.bwt.next_chars[c]] - 1

        div = 1./self.compression_occ

        while ep > sp and i >= 1:
            c = sample[i - 1]

            if sp % self.compression_occ == 0:
                rank_sp = self.bwt.tally[c][int(sp * div)]

            elif sp % self.compression_occ < half_compression or sp > int(n * div) * self.compression_occ:
                count = 0
                for up in range(sp, sp - (sp % self.compression_occ), -1):
                    if self.bwt.code[up] == c:
                        count += 1

                rank_sp = self.bwt.tally[c][int(sp * div)] + count

            else:
                count = 0
                for down in range(sp + 1, sp + (self.compression_occ - (sp % self.compression_occ)) + 1):
                    if self.bwt.code[down] == c:
                        count += 1

                rank_sp = self.bwt.tally[c][int(sp * div + 1)] - count

            if ep % self.compression_occ == 0:
                rank_ep = self.bwt.tally[c][int(ep * div)]

            elif ep % self.compression_occ < half_compression or ep > int(n * div) * self.compression_occ:
                count = 0
                for up in range(ep, ep - (ep % self.compression_occ), -1):
                    if self.bwt.code[up] == c:
                        count += 1

                rank_ep = self.bwt.tally[c][int(ep * div)] + count

            else:
                count = 0
                for down in range(ep + 1, ep + (self.compression_occ - (ep % self.compression_occ)) + 1):
                    if self.bwt.code[down] == c:
                        count += 1
                rank_ep = self.bwt.tally[c][int(ep * div + 1)] - count

            sp = self.bwt.f[c] + rank_sp - 1
            ep = self.bwt.f[c] + rank_ep - 1
            i = i - 1

        if sp >= ep:
            return []
        else:
            return [self.bwt.get_sa(i) for i in range(sp + 1, ep + 1)]

    def write(self, path_to_file: str = "index.data") -> None:
        """
        writing index_data to disk (index.data)
        """

        file = bz2.BZ2File(path_to_file, 'w')
        pickle.dump(self, file)
        file.close()

    @classmethod
    def read(cls, path_to_file: str = "index.data") -> 'FmIndex':
        """
        returns index that is stored in path
        """

        file = bz2.BZ2File(path_to_file, 'r')
        index = pickle.load(file)
        file.close()

        return index


if __name__ == "__main__":
    ref_genome = 'TTTAAAGGGNNNCAG'

    sample = 'TTTAAA'

    print("Reference genome: ", ref_genome)
    print("Sample: ", sample, "\n")
    print("len ", len(ref_genome))
    index = FmIndex(ref_genome)

    print(F'kmer match {sample}')
    match = index.query(sample)
    print(match, '\n')

    for m in match:
        assert (ref_genome[m:(m + len(sample))] == sample)

    print('perfect match')
    perfect_match = index.query(sample)

    print(perfect_match)

    code = index.bwt.get_bwt(ref_genome)

    print(code)
    print(str(index))

    index.write()

    index2 = FmIndex.read()

    print(F'kmer match {sample}')
    match = index2.query(sample)
    print(match, '\n')

    for m in match:
        assert (ref_genome[m:(m + len(sample))] == sample)

    print('perfect match')
    perfect_match = index2.query(sample)

    print(perfect_match)

    code = index2.bwt.get_bwt(ref_genome)
    print(code)
    print(str(index2))

    print("compression")

    index_compressed2 = FmIndex(ref_genome, compression_occ=9)

    match = index_compressed2.query(sample)
    print(match)


    index_compressed2.query(sample)

    print("TEST")
    index = FmIndex(ref_genome, compression_sa=1)
    start = time.perf_counter_ns()
    index.query("TTT")
    end = time.perf_counter_ns()

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/index_data/index.data")

    # The following strings are copied from the original genome


    ns = 10**-9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                                       "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                                       "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns*(end - start))

    print(getsizeof(index.bwt))