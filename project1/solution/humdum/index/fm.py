# LB, pre- 2020-10-09


# TODO:
# -tally matrix can be compressed
# -BWT can be compressed
# -suffix array can be compressed

# DONE:
# First column compressed to a dict of 6 entries

"""
Paolo Ferragina and Giovanni Manzini
An experimental study of a compressed index_data
InformationSciences, 135(1-2):13–28, June 2001
ISSN 00200255
doi:10.1016/S0020-0255(01)00098-6
URL: http://linkinghub.elsevier.com/retrieve/pii/S0020025501000986.
"""

# How the FM-Index works by Ben Langmead
# https://www.youtube.com/watch?v=kvVGj5V65io
# https://www.youtube.com/watch?v=UHGgpfxlaiE
# https://www.youtube.com/watch?v=DYyXgxYmSYU


from humdum.index import BurrowsWheeler
from typing import List
import bz2
import pickle

from humdum.utils import minidict


class FmIndex:
    """
    Implements the FM-index_data for substring lookup.

    Upon creation
        Burrows Wheeler (BW) object of size: 2 x len(reference genome)
        compressed version of the first column in the BW matrix (size: dict of 6 entries)
        tally/rank/occurrence matrix of size: 5 x len(reference_genome)
        helper dict to find the next lexicographical character in the context of DNA (size: dict of 5 entries)
    is stored.

    The FM-index_data works only for strings over the alphabet {A, C, G, T}.
    """

    def __string_checks(self, string):
        if ('$' in string):
            raise ValueError("The string should not contain '$'.")

        if (not set(string).issubset(set("ACGTN"))):
            raise ValueError("The string contains unexpected characters.")

        if (not string):
            raise ValueError("The string may not be empty.")

    def __init__(self, reference_genome: str):
        self.__string_checks(reference_genome)

        # helper dictionary to determine the next character (used in query method)
        self.next_chars = minidict({'$': 'A', 'A': 'C', 'C': 'G', 'G': 'N', 'N': 'T', 'T': None})

        # This creates a copy of the whole string?
        reference_genome += "$"

        # contains the Burrows Wheeler transformation of the reference genome
        # and the offsets of the suffix array
        # so far uncompressed => both of size len(ref_genome)
        self.bwt = BurrowsWheeler(reference_genome)

        # compressed first column of Burrows Wheeler matrix (e.g. cumulative frequencies of characters)
        self.F = self._shifts_F(reference_genome)

        # ranks of bases / occurrence matrix
        # Structure: A | C | G | N | T
        self.tally = self._build_tally(self.bwt.code, 1)

    def _shifts_F(self, reference_genome: str) -> dict:
        """
        Returns the shifts in the first column of the Burrows Wheeler matrix (compressed).
        Works only for DNA strings over the alphabet {A,C,G,T}.
        """

        # `None` key is added later
        shifts = {k: 0 for k in self.next_chars.keys()}

        for i in range(len(reference_genome)):
            shifts[reference_genome[i]] = shifts[reference_genome[i]] + 1

        count_a = shifts['$']
        count_c = count_a + shifts['A']
        count_g = count_c + shifts['C']
        count_n = count_g + shifts['G']
        count_t = count_n + shifts['N']

        shifts['$'] = 0
        shifts['A'] = count_a
        shifts['C'] = count_c
        shifts['G'] = count_g
        shifts['N'] = count_n
        shifts['T'] = count_t
        shifts[None] = len(reference_genome)

        return shifts

    def _build_tally(self, bw_transform: str, step: int = 1) -> dict:
        """
        Returns tallies, i.e. ranks/occurrence of characters.
        """

        s = [int(bw_transform[0] == '$')]
        a = [int(bw_transform[0] == 'A')]
        c = [int(bw_transform[0] == 'C')]
        g = [int(bw_transform[0] == 'G')]
        n = [int(bw_transform[0] == 'N')]
        t = [int(bw_transform[0] == 'T')]

        count_s = s[0]
        count_a = a[0]
        count_c = c[0]
        count_g = g[0]
        count_n = n[0]
        count_t = t[0]

        for (count, item) in enumerate(bw_transform[1:], start=1):
            count_s = count_s + (item == '$')
            count_a = count_a + (item == 'A')
            count_c = count_c + (item == 'C')
            count_g = count_g + (item == 'G')
            count_n = count_n + (item == 'N')
            count_t = count_t + (item == 'T')

            if not (count % step):
                s.append(count_s)
                a.append(count_a)
                c.append(count_c)
                g.append(count_g)
                n.append(count_n)
                t.append(count_t)

        return {'$': s, 'A': a, 'C': c, 'G': g, 'N': n, 'T': t}

    def __len__(self):
        """
        The length of the reference string.
        """
        return (len(self.bwt.code) - 1)

    def __str__(self):
        """
        Returns the original string.
        """
        return self.bwt.decode()

    def query(self, sample: str) -> List[int]:
        """
        Query index_data for a given sample.
        Returns a list of all positions in the reference genome (0-based).
        """

        self.__string_checks(sample)

        if (len(sample) >= len(self.bwt.code)):
            raise ValueError(
                F"Sample length may not exceed that of the reference genome ({len(self)})."
            )

        i = len(sample) - 1
        c = sample[i]
        sp = self.F[c]
        ep = self.F[self.next_chars[c]] - 1

        while (ep >= sp) and (i >= 1):
            c = sample[i - 1]
            sp = self.F[c] + self.tally[c][sp - 1]
            ep = self.F[c] + self.tally[c][ep] - 1
            i = i - 1

        if (sp > ep):
            return []
        else:
            return [self.bwt.sa[i] for i in range(sp, ep + 1)]

    def write(self, path: str = "", file_name: str = "index.data") -> None:
        """
        writing index_data to disk (index.data)
        """

        full_name = path + file_name

        file = bz2.BZ2File(full_name, 'w')
        pickle.dump(self, file)
        file.close()

    @classmethod
    def read(cls, path: str = "", file_name: str = "index.data") -> 'FmIndex':
        """
        returns index that is stored in path
        """

        full_name = path + file_name

        file = bz2.BZ2File(full_name, 'r')
        index = pickle.load(file)
        file.close()

        return index


if __name__ == "__main__":
    ref_genome = 'TAGAGAGATCGATCGACTGACTGACTCAGN'

    sample = 'CTCAGN'

    print("Reference genome: ", ref_genome)
    print("Sample: ", sample, "\n")

    index = FmIndex(ref_genome)

    print(F'kmer match {sample}')
    match = index.query(sample)
    print(match, '\n')

    for m in match:
        assert (ref_genome[m:(m + len(sample))] == sample)

    print('perfect match')
    perfect_match = index.query(sample)

    print(perfect_match)

    print(index.bwt.code)
    print(index.bwt.decode())

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

    print(index2.bwt.code)
    print(index2.bwt.decode())
