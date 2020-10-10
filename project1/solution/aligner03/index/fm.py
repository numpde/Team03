# LB, pre- 2020-10-09


# TODO:
# -tally matrix can be compressed
# -BWT can be compressed
# -suffix array can be compressed

# DONE:
# First column compressed to a dict of 6 entries


from aligner03.index import BurrowsWheeler
from typing import List


class FmIndex:
    """
    Can be used to create a FM-index.
    Upon creation
        Burrows Wheeler (BW) object (size: 2* len(reference genome))
        compressed version of the first column in the BW matrix (size: dict of 6 entries)
        tally/rank/occurrence matrix (size: 5*len(reference_genome)
        helper dict to find the next lexicographical character in the context of DNA (size: dict of 5 entries)
    is stored.
    The FM-index works only for strings that contains the characters ['A','C','G','T']
    """

    def __string_checks(self, string):
        if ('$' in string):
            raise ValueError("The string should not contain '$'.")

        if (not set(string).issubset(set("ACGTN"))):
            raise ValueError("The string contains unexpected characters.")

        if ('N' in string):
            raise NotImplementedError("The string contains the character 'N'.")

        if (not string):
            raise ValueError("The string may not be empty.")

    def __init__(self, reference_genome: str):
        self.__string_checks(reference_genome)

        # This creates a copy of the whole string?
        reference_genome += "$"

        # contains the Burrows Wheeler transformation of the reference genome
        # and the offsets of the suffix array
        # so far uncompressed => both of size len(ref_genome)
        self.bwt = BurrowsWheeler(reference_genome)

        # compressed first column of Burrows Wheeler matrix (e.g. cumulative frequencies of characters)
        self.F = self._shifts_F(reference_genome)

        # ranks of bases/ occurrence matrix
        # Structure: A| C | G | T
        self.tally = self._build_tally(self.bwt.code, 1)

        # helper dictionary to determine the next character (used in query method)
        self.next_chars = {'$': 'A', 'A': 'C', 'C': 'G', 'G': 'T', 'T': None}

    def _shifts_F(self, reference_genome: str) -> dict:
        """
        Returns the shifts in the first column of the Burrows Wheeler matrix (compressed).
        Works only for DNA strings over the alphabet {A,C,G,T}.
        """

        shifts = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0, None: 0}

        for i in range(len(reference_genome)):
            shifts[reference_genome[i]] = shifts[reference_genome[i]] + 1

        count_A = shifts['$']
        count_C = count_A + shifts['A']
        count_G = count_C + shifts['C']
        count_T = count_G + shifts['G']

        shifts['$'] = 0
        shifts['A'] = count_A
        shifts['C'] = count_C
        shifts['G'] = count_G
        shifts['T'] = count_T
        shifts[None] = len(reference_genome)

        return shifts

    def _build_tally(self, bw_transform: str, step: int = 1) -> dict:
        """
        Returns tallies, i.e. ranks/occurrence of characters.
        """

        S = [int(bw_transform[0] == '$')]
        A = [int(bw_transform[0] == 'A')]
        C = [int(bw_transform[0] == 'C')]
        G = [int(bw_transform[0] == 'G')]
        T = [int(bw_transform[0] == 'T')]

        count_S = S[0]
        count_A = A[0]
        count_C = C[0]
        count_G = G[0]
        count_T = T[0]

        for (count, item) in enumerate(bw_transform[1:], start=1):
            count_S = count_S + (item == '$')
            count_A = count_A + (item == 'A')
            count_C = count_C + (item == 'C')
            count_G = count_G + (item == 'G')
            count_T = count_T + (item == 'T')

            if not (count % step):
                S.append(count_S)
                A.append(count_A)
                C.append(count_C)
                G.append(count_G)
                T.append(count_T)

        return {'$': S, 'A': A, 'C': C, 'G': G, 'T': T}

    def __len__(self):
        return (len(self.bwt.code) - 1)

    def query(self, sample: str) -> List[int]:
        """
        Query index for a given sample.
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


if __name__ == "__main__":
    ref_genome = 'TAGAGAGATCGATCGACTGACTGACTCAG'

    sample = 'ACT'

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
