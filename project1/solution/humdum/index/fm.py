# LB, pre- 2020-10-09

import time
from humdum.index import BurrowsWheeler
from humdum.index.wt import WaveletTree
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
            raise ValueError("The string must not be empty.")

    def __init__(self, reference_genome: str, compression_occ: int = 32, compression_sa: int = 32, wavelet=True):

        if compression_occ < 1 or compression_sa < 1:
            raise ValueError("compression coefficients need to be strictly positive >=0")

        self.__string_checks(reference_genome)

        self.compression_occ = compression_occ

        self.n = len(reference_genome)

        reference_genome += "$"

        # contains the Suffix Array of the reference genome
        # can create the Burrows Wheeler Transformation on demand
        # size: len(ref_genome) * sizeof(int)

        if wavelet:
            self.bwt = WaveletTree(reference_genome, compression_sa=compression_sa)
        else:
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

    def query(self, sample: str) -> List[int]:
        """
        Query index_data for a given sample.
        Returns a list of all positions in the reference genome (0-based).
        """

        self.__string_checks(sample)

        if len(sample) > self.n:
            raise ValueError(
                F"Sample length may not exceed that of the reference genome ({len(self)})."
            )

        n = self.n

        i = len(sample) - 1
        c = sample[i]
        sp = self.bwt.f[c] - 1
        ep = self.bwt.f[self.bwt.next_chars[c]] - 1

        while ep > sp and i >= 1:
            c = sample[i - 1]

            rank_sp = self.bwt.rank(c, sp)
            rank_ep = self.bwt.rank(c, ep)

            sp = self.bwt.f[c] + rank_sp - 1
            ep = self.bwt.f[c] + rank_ep - 1

            i = i - 1

        if sp >= ep:
            return []
        else:
            return [self.bwt.get_sa(i) for i in range(sp + 1, ep + 1)]

    def write(self, path_to_file):
        """
        Write index to file.

        Returns self.
        """

        with bz2.BZ2File(str(path_to_file), mode='w') as fd:
            pickle.dump(self, fd)

        return self

    @classmethod
    def read(cls, path_to_file) -> 'FmIndex':
        """
        Recover the index from file.
        """

        with bz2.BZ2File(str(path_to_file), 'r') as fd:
            index = pickle.load(fd)
            import humdum
            assert isinstance(index, FmIndex) or isinstance(index, humdum.index.FmIndex)
            return index

    @classmethod
    def read_or_make(cls, *, path_to_genome, path_to_index=None):
        """
        Create an index for the genome and write to file.
        Attempt to read from file instead if it already exists.
        Default `path_to_index` appends the suffix ".index".
        Returns the index.

        RA, 2020-10-23
        """

        from pathlib import Path
        DEFAULT_SUFFIX = ".index"
        path_to_genome = Path(path_to_genome)
        path_to_index = Path(path_to_index or (str(path_to_genome) + DEFAULT_SUFFIX))

        assert path_to_genome.is_file()

        if path_to_index.is_file():
            return cls.read(path_to_index)
        else:
            from humdum.io import from_fasta
            from humdum.utils import unlist1
            return cls(unlist1(list(from_fasta(path_to_genome))).seq).write(path_to_index)

    def query_hist(self, sample: str, hist: List[float]) -> List[float]:
        """
        Returns the time it takes to query the sample
        Result of the sample query is not returned (use method query(str))
        """
        start = time.perf_counter_ns()
        self.query(sample)
        end = time.perf_counter_ns()
        hist.append(end - start)
        return hist


def plot_hist(hist: List[float]):
    """
    Plots a histogram around the mean +- 2 standard deviations
    """
    import numpy as np
    import matplotlib.pyplot as plt

    hist = np.array(hist)

    mu, sigma = hist.mean(), hist.var()

    # the histogram of the data
    n, bins, patches = plt.hist(hist, bins=int(np.sqrt(len(hist))),
                                range=(mu - 2 * sigma ** 0.5, mu + 2 * sigma ** 0.5), density=True)

    plt.xlabel('Time [ns]')
    plt.ylabel('Frequency')
    plt.title(r'$\mathrm{Histogram\ of\ query time:}\ \mu= %1.f \ ns,\ \sigma= %1.f \ ns$' % (mu, sigma))
    plt.grid(True)

    plt.show()


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

    index.write("index.data")

    index2 = FmIndex.read("index.data")

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
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    hist1 = []
    hist2 = []
    hist3 = []
    hist4 = []

    for i in range(10000):
        hist1 = index.query_hist("AAAAGAATGCA", hist1)

    for i in range(10000):
        hist2 = index.query_hist("CGACACCACCAAGGCCACCCACCTGCCT", hist2)

    for i in range(10000):
        hist3 = index.query_hist("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG", hist3)

    for i in range(10000):
        hist4 = index.query_hist("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                                 "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA", hist4)

    plot_hist(hist1)
    plot_hist(hist2)
    plot_hist(hist3)
    plot_hist(hist4)

    hist_test = []
    for i in range(10000):
        hist_test = index.query_hist("AAAAGAATGCA", hist_test)
        hist_test = index.query_hist("CGACACCACCAAGGCCACCCACCTGCCT", hist_test)
        hist_test = index.query_hist("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG", hist_test)
        hist_test = index.query_hist("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                                     "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA", hist_test)

    plot_hist((hist_test))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.sa16_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.sa16_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet8_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet16_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet32_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet8_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read wavelet1")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet1_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))

    print("read")
    index = FmIndex.read("../../tests/data_for_tests/data/genome.chr22.fa.gz.wavelet32_index")

    # The following strings are copied from the original genome

    ns = 10 ** -9

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("CGACACCACCAAGGCCACCCACCTGCCT"))
    start = time.perf_counter_ns()
    index.query("CGACACCACCAAGGCCACCCACCTGCCT")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG"))
    start = time.perf_counter_ns()
    index.query("GGCATTTACAACTAAAACATTGAATTCAGATTCATTTTCAGGTAATGATATAATCATGTG")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    print(len("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
              "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA"))
    start = time.perf_counter_ns()
    index.query("AAAAGAATGCATTTCTGTATTTTTTGAAACCTTTTCTTTTGAAAACATAGTAATACATTT"
                "CTACTCTAAAATAGAACTTAGCCTAAATACTTTCAAAACCTTTAGAATTTGGAAAAGAAA")
    end = time.perf_counter_ns()
    print("time: ", ns * (end - start))

    # print(getsizeof(index.bwt))
