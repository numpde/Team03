import sys
from typing import List, Tuple
from BitVector import BitVector
from humdum.utils import minidict

from objsize import get_deep_size


class WaveletTree:

    def __init__(self, reference_genome: str ):
        # Bitvector of nodes; size = 2.5 * len(ref_genome)
        self.bits = self.create_bit_vecs(reference_genome)

        # Structure [Parent, Me, Left_child, Right_Child]
        self.meta = [[None, 0, 1, 2], [0, 1, 'N', 'A'], [0, 2, 3, 4], [2, 3, 'C', 'G'], [2, 4, 'T', '$']]

        # 0 = left, 1 = right
        self.codes = minidict({'N': [0, 0], 'A': [0, 1], 'C': [1, 0, 1], 'G': [1, 0, 1], 'T': [1, 1, 0], '$': [1, 1, 1]})

    def create_bit_vecs(self, reference_genome:str) -> List[BitVector]:

        #Placeholders
        bit_vec0 = BitVector(bitlist=[1, 0, 0, 1, 1, 0, 0, 1])

        bit_vec1 = BitVector(bitlist=[1, 0, 0, 1])
        bit_vec2 = BitVector(bitlist=[1, 1, 0, 0])

        bit_vec3 = BitVector(bitlist=[1, 1])
        bit_vec4 = BitVector(bitlist=[1, 0])

        return [bit_vec0, bit_vec1, bit_vec2, bit_vec3, bit_vec4]

    def get_root(self) -> int:
        return 0

    def get_left_child(self, index: int) -> int:
        return self.meta[index][2]

    def get_right_child(self, index: int) -> int:
        return self.meta[index][3]

    def get_parent(self, index: int) -> int:
        return self.meta[index][0]

    def rank(self, char: str, index: int) -> int:
        return 0


if __name__ == "__main__":
    ref_genome = 'ACAAACGTACGATCGACTACTACA'
    sample = 'CTCAGN'

    print("Reference genome: ", ref_genome)
    print("Sample: ", sample, "\n")

    tree = WaveletTree(ref_genome)

    print(tree.meta[0])

    print(tree.get_left_child(tree.get_root()))
    print(tree.get_right_child(tree.get_root()))
    print(tree.get_parent(4))