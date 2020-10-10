# LB, pre- 2020-10-09
# RA, 2020-10-09

from typing import List, Tuple


class BurrowsWheeler:
    """
    Generates the Burrows Wheeler transform (BWT) and
    the corresponding offsets in the suffix array (OSA).
    Upon creation of the object the BWT and OSA
    are created and stored (Both of size len(reference_genome)).
    A method to decode a BWT is provided.
    Works with strings of any shape.
    The encoding can be done via suffix arrays (default)
    or via rotation of the string (specified at initialisation).
    """

    def __init__(self, reference_genome: str, strategy: str = 'S'):
        # choose encoding strategy to build the Burrows-Wheeler Transform
        # R = via rotation of the string
        # S = via suffix array (preferred method)

        if strategy not in ['S', 'R']:
            raise ValueError('strategy needs to be S (suffix array) or R (rotation)')

        if len(reference_genome) < 1:
            raise ValueError('please provide a non empty string for the reference genome')


        self.encode_strategy = strategy

        reference_genome = reference_genome + '$'

        # get Burrows Wheeler transformation and corresponding offsets in the suffix array
        (self.code, self.sa) = self.bwt_encode(reference_genome)

    def bwt_encode(self, reference_genome: str):
        '''
        Returns the Burrows Wheeler transformation, given a reference genome,
        '''

        if (self.encode_strategy == 'R'):
            return self.bwt_encode_rot(reference_genome)
        else:
            return self.bwt_encode_suffix(reference_genome)

    def bwt_encode_rot(self, reference_genome: str) -> str:
        '''
        Returns the Burrows-Wheeler Transform
        constructed by the rotation of the reference_genome, given a reference_genome
        '''

        n = len(reference_genome)
        bw_transform = sorted([reference_genome[i:n] + reference_genome[0:i] for i in range(n)])
        L = ''.join(last[-1] for last in bw_transform)

        return L

    def bwt_encode_suffix(self, reference_genome: str) -> Tuple[str, List[int]]:
        '''
        Returns the Burrows-Wheeler Transform (suffix, offset)
        constructed by the suffix array of the reference_genome
        '''

        suffix_array_offsets = self.suffix_array_offsets(reference_genome)
        bw_transform = []

        for w in suffix_array_offsets:
            if w == 0:
                bw_transform.append('$')
            else:
                bw_transform.append(reference_genome[w - 1])

        return (''.join(bw_transform), suffix_array_offsets)

    def suffix_array(self, reference_genome: str) -> List[Tuple[str, int]]:
        '''
        Returns the sorted suffix array (suffix, offset), given a reference_genome,
        '''

        n = len(reference_genome)
        suffix_array = sorted([(reference_genome[i:], i) for i in range(n)])

        return suffix_array

    def suffix_array_offsets(self, reference_genome):
        '''
        Returns the offsets of the sorted suffixes in reference genome
        Given a reference_genome
        '''

        suffix_array = self.suffix_array(reference_genome)
        offsets = [i[1] for i in suffix_array]

        return (offsets)

    def decode(self, bw_transform: str = None) -> str:
        '''
        Returns original string
        '''

        bw_transform = bw_transform or self.code

        matrix = [''] * len(bw_transform)

        # add Burrows Wheeler transform as first column and sort
        for i in range(len(bw_transform)):
            matrix = sorted(bw_transform[i] + matrix[i] for i in range(len(bw_transform)))

        s = ''
        for row in matrix:
            if row.endswith("$"):
                s = row

        return s[:-1]
