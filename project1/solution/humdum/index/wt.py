import sys
import time
from typing import List, Tuple
from BitVector import BitVector
from humdum.utils import minidict

from objsize import get_deep_size


class WaveletTree:
    """
    Generates a Suffix Array (SA), first column of Suffix Matrix (f)
    and a Wavelet Tree for the Burrows Wheeler transformation
    Upon creation of the object:
    - SA is created, compressed and stored
        len(reference_genome) / compression_sa * sizeof(int).
    - f is created and stored
        6 * sizeof(dict(char:int))
    - wavelet tree is created and stored
        2.5 len(reference_genome)
    - helper data structures
        ~ len(reference_genome)
    Works for strings over the alphabet {A, C, G, N, T}.
    The compression for the suffix array (compression_sa) can be select upon creation of the object.
    Different algorithms are provided and can be selected upon creation of the object:
    - KaerkkaeinenSanders:
        Runtime: O(n)
        Space: O(n)
    - ManberMyers:
        Runtime: O(nlog(n) )
        Space: O(n)
    - Simple:
        Runtime: O(n^2 log(n) )
        Space: O(n^2)
    """

    def __init__(self, reference_genome: str, strategy: str = 'KaerkkaeinenSanders',
                 compression_occ: int = 8, compression_sa: int = 8):

        if strategy not in ['KaerkkaeinenSanders', 'ManberMyers', 'Simple']:
            raise ValueError('strategy needs to be KaerkkaeinenSanders, ManberMyers or Simple ')

        if len(reference_genome) < 1:
            raise ValueError('please provide a non empty string for the reference genome')

        if compression_occ < 1 or compression_sa < 1:
            raise ValueError("compression coefficients need to be strictly positive >=0")

        if reference_genome[-1] != '$':
            reference_genome = reference_genome + '$'

        self.compression_occ = compression_occ
        self.compression_sa = compression_sa

        # helper dictionary to determine the next character (used in query method)
        self.next_chars = minidict({'$': 'A', 'A': 'C', 'C': 'G', 'G': 'N', 'N': 'T', 'T': None})

        # get Burrows Wheeler transformation and corresponding offsets in the suffix array
        self.sa, self.bitvector, bwt = self.suffix_array(reference_genome, strategy, compression_sa)

        # compressed first column of Burrows Wheeler matrix (e.g. cumulative frequencies of characters)
        self.f = self._shifts_f(reference_genome)

        # Bitvector of nodes; size = 2.5 * len(bwt)

        self.bits = self.create_bit_vecs(bwt)

        # Structure [Parent, Me, Left_child, Right_Child]
        self.meta = [[None, 0, 1, 2], [0, 1, 'N', 'A'], [0, 2, 3, 4], [2, 3, 'C', 'G'], [2, 4, 'T', '$'],
                     [1, 'N', None, None], [1, 'A', None, None], [1, 'C', None, None], [1, 'G', None, None],
                     [1, 'T', None, None], [1, '$', None, None]]

        # 0 = left, 1 = right
        self.codes = minidict({'N': [0, 0], 'A': [0, 1], 'C': [1, 0, 0],
                               'G': [1, 0, 1], 'T': [1, 1, 0], '$': [1, 1, 1]})

        self.n = len(reference_genome)

    def _shifts_f(self, reference_genome: str) -> dict:
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

    def suffix_array(self, reference_genome: str, strategy: str,
                     compression: int = 1) -> Tuple[List[int], BitVector, str]:
        """
        Returns the Burrows-Wheeler Transform (suffix, offset)
        constructed by the suffix array of the reference_genome
        """

        suffix_array = []
        if strategy == 'Simple':
            suffix_array = self.suffix_array_simple(reference_genome)
        elif strategy == 'ManberMyers':
            suffix_array = self.suffix_array_manbermyers(reference_genome)
        else:
            suffix_array = self.suffix_array_kaerkkaeinensanders(reference_genome, len(reference_genome), 6)

        code = self.get_bwt(reference_genome, suffix_array)
        return ([num for num in suffix_array if num % compression == 0],
                BitVector(bitlist=[1 if num % compression == 0 else 0 for num in suffix_array]), code)

    def suffix_array_kaerkkaeinensanders(self, reference_genome, n: int, k: int) -> List[int]:

        def to_int(string: str) -> List[int]:
            str_to_int = {'$': 1, 'A': 2, 'C': 3, 'G': 4, 'N': 5, 'T': 6}
            return [str_to_int[char] for char in string]

        def leq2(a1: int, a2: int, b1: int, b2: int) -> bool:
            return a1 < b1 or a1 == b1 and a2 <= b2

        def leq3(a1: int, a2: int, a3: int, b1: int, b2: int, b3: int) -> bool:
            return a1 < b1 or a1 == b1 and leq2(a2, a3, b2, b3)

        def radix_pass(a: List[int], r: List[int], len_a: int, len_k: int) -> List[int]:
            c = [0] * (len_k + 1)
            b = [0] * len_a

            for i in range(len_a):
                c[r[a[i]]] += 1
            sum = 0
            for i in range(len_k + 1):
                t = c[i]
                c[i] = sum
                sum += t
            for i in range(len_a):
                b[c[r[a[i]]]] = a[i]
                c[r[a[i]]] += 1
            return b

        sa = [0] * n
        s = reference_genome
        if type(s[0]) is not int:
            s = to_int(reference_genome)
            s += [0, 0, 0]

        n0 = int((n + 2) / 3)
        n1 = int((n + 1) / 3)
        n2 = int(n / 3)
        n02 = n0 + n2

        s12 = [0] * (n02 + 3)
        s0 = [0] * n0

        j = 0
        for i in range(n + n0 - n1):
            if i % 3 != 0:
                s12[j] = i
                j += 1

        sa12 = radix_pass(s12, s[2:], n02, k) + [0, 0, 0]
        s12 = radix_pass(sa12, s[1:], n02, k) + [0, 0, 0]
        sa12 = radix_pass(s12, s, n02, k) + [0, 0, 0]

        name = 0
        c0 = -1
        c1 = -1
        c2 = -1
        for i in range(n02):
            if s[sa12[i]] != c0 or s[sa12[i] + 1] != c1 or s[sa12[i] + 2] != c2:
                name += 1
                c0 = s[sa12[i]]
                c1 = s[sa12[i] + 1]
                c2 = s[sa12[i] + 2]
            if sa12[i] % 3 == 1:
                s12[int(sa12[i] / 3)] = name
            else:
                s12[int(sa12[i] / 3) + n0] = name

        if name < n02:
            sa12 = self.suffix_array_kaerkkaeinensanders(s12, n02, name)
            for i in range(n02):
                s12[sa12[i]] = i + 1
        else:
            for i in range(n02):
                sa12[s12[i] - 1] = i

        j = 0
        for i in range(n02):
            if sa12[i] < n0:
                s0[j] = 3 * sa12[i]
                j += 1

        sa0 = radix_pass(s0, s, n0, k)

        p = 0
        t = n0 - n1
        k = 0

        while k < n:

            i = sa12[t] * 3 + 1 if sa12[t] < n0 else (sa12[t] - n0) * 3 + 2
            j = sa0[p]
            if leq2(s[i], s12[sa12[t] + n0], s[j], s12[int(j / 3)]) if sa12[t] < n0 \
                    else leq3(s[i], s[i + 1], s12[sa12[t] - n0 + 1], s[j], s[j + 1], s12[int(j / 3) + n0]):
                sa[k] = i
                t += 1
                if t == n02:
                    k += 1
                    while p < n0:
                        sa[k] = sa0[p]
                        p += 1
                        k += 1

            else:
                sa[k] = j
                p += 1
                if p == n0:
                    k += 1
                    while t < n02:
                        sa[k] = sa12[t] * 3 + 1 if sa12[t] < n0 else (sa12[t] - n0) * 3 + 2
                        k += 1
                        t += 1
            k += 1

        return sa

    def suffix_array_manbermyers(self, reference_genome: str) -> List[int]:

        def sort_chars(reference_genome: str) -> List[int]:
            n = len(reference_genome)
            order = [0] * n
            count = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'N': 0, 'T': 0}

            for i in range(n):
                count[reference_genome[i]] = count[reference_genome[i]] + 1
            keys = list(count.keys())
            for (i, j) in enumerate(count):
                if j == '$':
                    continue
                count[j] = count[j] + count[keys[i - 1]]

            for i in range(n - 1, -1, -1):
                c = reference_genome[i]
                count[c] = count[c] - 1
                order[count[c]] = i

            return order

        def compute_classes(reference_genome: str, order: List[int]) -> List[int]:
            n = len(reference_genome)
            classes = [0] * n
            classes[order[0]] = 0

            for i in range(1, n):
                if reference_genome[order[i]] != reference_genome[order[i - 1]]:
                    classes[order[i]] = classes[order[i - 1]] + 1
                else:
                    classes[order[i]] = classes[order[i - 1]]
            return classes

        def sort_doubled(reference_genome: str, step: int, order: List[int], classes: List[int]) -> List[int]:
            n = len(reference_genome)
            count = [0] * n
            new_order = [0] * n

            for i in range(0, n):
                count[classes[i]] = count[classes[i]] + 1
            for j in range(1, n):
                count[j] = count[j] + count[j - 1]
            for i in range(n - 1, -1, -1):
                start = (order[i] - step + n) % n
                cl = classes[start]
                count[cl] = count[cl] - 1
                new_order[count[cl]] = start

            return new_order

        def updated_classes(order: List[int], classes: List[int], step: int) -> List[int]:
            n = len(order)
            new_classes = [0] * n
            new_classes[order[0]] = 0
            for i in range(1, n):
                cur = order[i]
                prev = order[i - 1]
                mid = cur + step
                mid_prev = (prev + step) % n

                if classes[cur] != classes[prev] or classes[mid] != classes[mid_prev]:
                    new_classes[cur] = new_classes[prev] + 1
                else:
                    new_classes[cur] = new_classes[prev]

            return new_classes

        order = sort_chars(reference_genome)
        classes = compute_classes(reference_genome, order)

        step = 1
        n = len(reference_genome)
        while step < n:
            order = sort_doubled(reference_genome, step, order, classes)
            classes = updated_classes(order, classes, step)
            step *= 2

        return order

    def suffix_array_simple(self, reference_genome):
        """
        Returns the offsets of the sorted suffixes in reference genome
        Given a reference_genome
        """

        n = len(reference_genome)

        suffix_array = sorted([(reference_genome[i:], i) for i in range(n)])

        offsets = [i[1] for i in suffix_array]

        return offsets

    def get_bwt(self, reference_genome: str, suffix_array: List[int] = None) -> str:

        suffix_array = suffix_array or self.sa

        bw_transform = []

        for w in suffix_array:
            if w == 0:
                bw_transform.append('$')
            else:
                bw_transform.append(reference_genome[w - 1])

        return ''.join(bw_transform)

    def __len__(self):
        return self.n - 1

    def __str__(self):
        """
        return original string
        """
        half_compression = self.compression_occ * 0.5
        n = self.n - 1

        next_char = self.access(0)
        next_row = 0
        original = next_char

        rank = 0

        for i in range(n - 1):
            rank = self.rank(next_char, next_row)

            skip = self.f[next_char]
            next_row = rank + skip - 1

            next_char = self.access(next_row)
            original = next_char + original

        return original

    def __sizeof__(self):

        print("sizes:")
        print("compression_occ:\t ", get_deep_size(self.compression_occ))
        print("compression_sa:\t\t ", get_deep_size(self.compression_sa))
        print("next_chars\t\t\t ", get_deep_size(self.next_chars))
        print("SA:\t\t\t\t\t ", get_deep_size(self.sa))
        print("F:\t\t\t\t\t ", get_deep_size(self.f))
        print("bitvec:\t\t\t\t ", get_deep_size(self.bitvector))
        print("bits:\t\t\t\t ", get_deep_size(self.bits))
        print("meta:\t\t\t\t ", get_deep_size(self.meta))
        print("codes:\t\t\t\t ", get_deep_size(self.codes))
        print("n:\t\t\t\t\t ", get_deep_size(self.n))

        total = get_deep_size(self.compression_occ) + get_deep_size(self.compression_sa) + \
            get_deep_size(self.next_chars) + get_deep_size(self.sa) + \
            get_deep_size(self.f) + get_deep_size(self.bitvector) + get_deep_size(self.bits) + \
            get_deep_size(self.meta) + get_deep_size(self.codes) + get_deep_size(self.n)

        print("Total:\t\t\t\t ", total)

        return total

    def get_sa(self, index: int) -> int:
        """
        Return entry in Suffix Array at position index
        """

        if self.bitvector[index] == 1:
            return self.sa[self.bitvector.rank_of_bit_set_at_index(index) - 1]
        else:

            half_compression = self.compression_occ * 0.5
            n = self.n - 1

            next_char = self.access(index)
            next_row = index

            rank = 0
            counter = 0
            while self.bitvector[next_row] != 1:
                rank = self.rank(next_char, next_row)

                skip = self.f[next_char]
                next_row = rank + skip - 1
                next_char = self.access(next_row)

                counter += 1

            return self.sa[self.bitvector.rank_of_bit_set_at_index(next_row) - 1] + counter

    def create_bit_vecs(self, lbwt: str) -> List[BitVector]:

        bit_vec0 = BitVector(bitlist=[0 if char == 'N' or char == 'A' else 1 for char in lbwt])

        rbwt = [char for char in lbwt if char != 'N' and char != 'A']
        lbwt = [char for char in lbwt if char == 'N' or char == 'A']

        bit_vec1 = BitVector(bitlist=[0 if char == 'N' else 1 for char in lbwt])
        bit_vec2 = BitVector(bitlist=[0 if char == 'C' or char == 'G' else 1 for char in rbwt])

        lbwt = [char for char in rbwt if char == 'C' or char == 'G']
        rbwt = [char for char in rbwt if char == 'T' or char == '$']

        bit_vec3 = BitVector(bitlist=[0 if char == 'C' else 1 for char in lbwt])
        bit_vec4 = BitVector(bitlist=[0 if char == 'T' else 1 for char in rbwt])

        return [bit_vec0, bit_vec1, bit_vec2, bit_vec3, bit_vec4]

    def get_root(self) -> int:
        return 0

    def get_left_child(self, node: int) -> int:
        return self.meta[node][2]

    def get_right_child(self, node: int) -> int:
        return self.meta[node][3]

    def get_parent(self, node: int) -> int:
        return self.meta[node][0]

    def rank(self, char: str, index: int):

        codes = self.codes[char]
        curr_node = 0
        curr_index = index
        rank = 0

        for code in codes:

            rank = self.rank_bit(curr_index, curr_node)

            if self.bits[curr_node][curr_index] != code:
                rank = curr_index - rank + 1

            if rank == 0:
                return 0

            curr_index = rank - 1

            if code == 0:
                curr_node = self.get_left_child(curr_node)
            else:
                curr_node = self.get_right_child(curr_node)

        return rank

    def rank_bit(self, index: int, node: int = 0) -> int:
        if self.bits[node][index] == 1:
            return self.bits[node].rank_of_bit_set_at_index(index)
        else:
            next_set_index = self.bits[node].next_set_bit(index)
            if next_set_index == -1:
                i = index
                while self.bits[node][i] != 1 and i >= 0:
                    i -= 1
                if i >= 0:
                    return i + 1 - self.bits[node].rank_of_bit_set_at_index(i) + (index - i)
                else:
                    return index + 1
            else:
                rank0 = next_set_index + 1 - self.bits[node].rank_of_bit_set_at_index(next_set_index)
            return rank0 - (next_set_index - index - 1)

    def access(self, index: int, node: int = 0) -> str:

        curr_node = node

        # leaf node
        if self.meta[curr_node][2] is None:
            return self.meta[curr_node][1]

        # parent of a leaf
        if type(self.meta[curr_node][2]) is str:
            return self.meta[curr_node][2]

        bit = self.bits[curr_node][index]
        curr_index = index

        while type(self.meta[curr_node][2]) is not str:

            par_node = curr_node
            if bit == 0:
                curr_node = self.get_left_child(curr_node)
            else:
                curr_node = self.get_right_child(curr_node)

            curr_index = self.rank_bit(curr_index, par_node) - 1
            bit = self.bits[curr_node][curr_index]

        return self.meta[curr_node][2+bit]


if __name__ == "__main__":

    ref_genome = 'AGCTA'

    sample = 'CTCAGN'

    print("Reference genome: ", ref_genome)
    print("Sample: ", sample, "\n")

    bwt = WaveletTree(ref_genome, strategy='Simple')

    print(bwt.sa)

    print(bwt.get_bwt(ref_genome))

    bwt = WaveletTree(ref_genome, strategy='ManberMyers', compression_sa=1)

    print(str(bwt))

    print(bwt.get_sa(1))

    print(bwt.sa)
    print(bwt.get_bwt(ref_genome))

    sys.getsizeof(bwt)

    bwt = WaveletTree(ref_genome, strategy='KaerkkaeinenSanders', compression_occ=10, compression_sa=10)

    print(bwt.sa)
    print(bwt.get_bwt(ref_genome))

    print(str(bwt.bitvector))

    print(str(bwt))
