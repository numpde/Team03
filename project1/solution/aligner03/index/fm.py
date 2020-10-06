# Can be used to create a FM-index. Upon creation a Burrows Wheeler (BW) object (size: 2* len(reference genome)), a compressed version of the first column in the BW matrix (size: dict of 6 entries), tally/rank matrix (size: 5*len(reference_genome), the kmer constant (size: 1 int) and a helper dict to find the next lexicographical character in the context of DNA (size: dict of 5 entries) is stored. The FM-index works only for strings that contains the characters ['A','C','G','T']


# TODO:
# -tally matrix can be compressed
# -BWT can be compressed
# -suffix array can be compressed

# DONE:
# First column compressed to a dict of 6 entries


from bw import Burrows_Wheeler


class FM_Index(object):

    def __init__(self, reference_genome, k=3):

        if (k < 1 or k > len(reference_genome)):
            raise IndexError(
                ' k has to be greater than 0 and smaller than the length of the reference genome. Make sure reference genome is not empty')

        # make sure the reference_genome terminates with a '$'
        if (reference_genome[-1] != '$'):
            reference_genome = reference_genome + "$"

        # contains the Burrows Wheeler transformation of the reference genome and the offsets of the suffix array (so far uncompressed => both of size len(ref_genome)
        self.bwt = Burrows_Wheeler(reference_genome)

        # compressed first column of Burrows Wheeler matrix (e.g. cumulative frequencies of characters)
        self.F = self.shifts_F(reference_genome)

        # number of characters used (kmer) to find a match in the index
        self.k = k

        # ranks of bases
        # Structure: A| C | G | T
        self.tally = self.build_tally(self.bwt.code, 1)

        # helper dictionary to determine the next character (used in query method)
        self.next_chars = {'$': 'A', 'A': 'C', 'C': 'G', 'G': 'T', 'T': 'total'}

    # Returns the shifts in the first column of the Burrows Wheeler matrix (compressed)
    # Works only for DNA strings (A,C,G,T)
    def shifts_F(self, reference_genome):

        shifts = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'total': 0}

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
        shifts['total'] = len(reference_genome)

        return shifts

    # Returns tally matrix (e.g. ranks of characters)
    def build_tally(self, bw_transform, step=1):

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

        for count, item in enumerate(bw_transform[1:], start=1):

            count_S = count_S + int(item == '$')
            count_A = count_A + int(item == 'A')
            count_C = count_C + int(item == 'C')
            count_G = count_G + int(item == 'G')
            count_T = count_T + int(item == 'T')

            if count % step == 0:
                S.append(count_S)
                A.append(count_A)
                C.append(count_C)
                G.append(count_G)
                T.append(count_T)

        return {'$': S, 'A': A, 'C': C, 'G': G, 'T': T}

    # query index for a given sample
    # this searches ONLY FOR A MATCH OF THE FIRST K characters of the sample
    # sample should be a string
    # return a list of all positions in the reference genome
    def query(self, sample, k=None):

        k = k or self.k

        if (len(sample) < k):
            raise ValueError('sample length needs to be greater than k = ' + str(self.k))

        if (len(sample) > len(self.bwt.code) - 1):
            raise ValueError(
                'sample length needs to be smaller than reference genome <= ' + str(len(self.bwt.code) - 1))

        kmer = sample[:k]

        i = len(kmer) - 1
        c = kmer[i]
        sp = self.F[c]
        ep = self.F[self.next_chars[c]] - 1

        while sp <= ep and i >= 1:
            c = sample[i - 1]
            sp = self.F[c] + self.tally[c][sp - 1]
            ep = self.F[c] + self.tally[c][ep] - 1
            i = i - 1

        if (sp > ep):
            return []
        else:
            return ([self.bwt.sa[i] for i in range(sp, ep + 1)])


# searching for a match of the first k letters of the sample in the reference genome
# sample should be string
def query_index(sample, index, k=None):
    return index.query(sample, k)


# searching for a perfect match of the sample in the reference genome
# sample should be string
def perfect_match(sample, index):
    return index.query(sample, len(sample))


if __name__ == "__main__":
    ref_genome = 'TAGAGAGATCGATCGACTGACTGACTCAG$'

    sample = 'ACTACTGTCA'

    print("reference genome: ", ref_genome)
    print("sample: ", sample, "\n")

    k = 3
    index = FM_Index(ref_genome, k)

    print('kmer match  (', sample[:k], ')')
    match = query_index(sample, index)
    print(match, '\n')

    print('perfect match')
    perfect_match = perfect_match(sample, index)

    print(perfect_match)
