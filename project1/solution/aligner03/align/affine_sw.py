# HK, pre- 2020-10-09
# RA, 2020-10-09

import itertools
import numpy as np
from typing import List, Tuple
import re

from frozendict import frozendict

default_mutation_costs = frozendict({
    # Deletion
    'D': -2,
    # Insertion
    'I': -2,
    # Mutation
    'X': -10,
    # Match
    '=': 10,
})


class Alignment:
    def __init__(self):
        self.cigar_string = ''
        self.start_pos = None  # 1-based start position of query in ref
        self.start_coord = None  # (i,j) coordinates of start_pos, 1-based
        self.end_coord = None  # (i,j) coordinates of end_pos, 1-based
        self.score = None

    def __repr__(self):
        return F"Aligns at {self.start_pos} with CIGAR = {self.cigar_string} and score = {self.score}"

    def prepend_to_cigar_string(self, symbol: str, count: int = 1):
        """
        Prepends a symbol to the cigar string.
        symbol: symbol to prepend
        count: count of symbol to prepend
        """
        if self.cigar_string and symbol == re.findall(r"[XIDS=]", self.cigar_string)[0]:
            current_count = re.findall(r"[0-9]+", self.cigar_string)[0]
            self.cigar_string = str(int(current_count) + count) + self.cigar_string[len(current_count):]
        else:
            self.cigar_string = f'{count}{symbol}' + self.cigar_string

    def count_total(self) -> int:
        """
        Returns the total count of matches, deletions, mutations and insertions. I.e. it returns the lengths of the
        expanded cigar string.
        Example: returns 2+1+3 = 6 for 2=1I3=
        """
        total = 0
        for count in re.findall(r"[0-9]+", self.cigar_string):
            total += int(count)
        return total

    def matching_subsegments(self) -> List[Tuple]:
        """
        Returns a list of all matching 1-based subsegments of the query.
        Example:
        AAAAC
        AAABC
        returns: [(1,3),(5,5)]
        """
        cigar = self.cigar_string
        idx = self.start_pos + 1  # +1 since start_pos is 0-based
        segment = []
        for (n, a) in re.findall(r"([0-9]+)([XIDS=])", cigar):
            n = int(n)
            if a == '=':
                segment.append((idx + 1, idx + n))
            idx += n

        return segment

    def visualize(self, *, ref: str, query: str):
        """
        Returns the reference and the query together with the expanded cigar string for viusalsation.
        """
        cigar = self.cigar_string
        total = self.count_total()
        end_pos = self.start_pos + total  # 0-based

        query = query[self.start_coord[1]:self.end_coord[1]]
        if self.start_pos:
            cigar = str(self.start_pos + 1) + 'S' + cigar
        if (len(ref) - 1 - end_pos) > 0:
            cigar += str(len(ref) - end_pos) + 'S'

        i = j = 0
        x = y = z = ""
        for (n, a) in re.findall(r"([0-9]+)([XIDS=])", cigar):
            n = int(n)
            z += a * n
            if (a in '=XDS'):
                x += ref[i:(i + n)]
                i += n
            if (a in '=XI'):
                y += query[j:(j + n)]
                j += n
            if (a == 'I'):
                x += "-" * n
            if (a == 'S'):
                y += " " * n
            if (a == 'D'):
                y += "-" * n
        return (x, y, z)


class SmithWaterman:
    """
    Example:
        ref = 'ATGGCCTC'
        query = 'ACGGCTC'
        aligner = SmithWaterman()
        for alignment in aligner(query=query, ref=ref):
            print(alignment.cigar_string)

        1=1X2=1I3=
    """

    def __init__(self, mutation_costs=default_mutation_costs):
        self.match_score = mutation_costs['=']
        self.mismatch_cost = mutation_costs['X']
        self.gap_cost = mutation_costs['D']
        self.insertion_cost = mutation_costs['I']
        self.opening_cost = -2
        self.extention_cost = -2

    def _compute_gap_costs(self, H, pos: Tuple, gap_type: str) -> Tuple:
        """
        Computes affine costs for gap extention.
        #TODO
        """
        i, j = pos
        if gap_type == 'D':
            costs = [a + c for a, c in
                     zip(H[:i, j][::-1],
                         range(self.opening_cost, self.opening_cost + len(H[:i, j]) * (self.extention_cost),
                               self.extention_cost))]
        elif gap_type == 'I':
            costs = [a + c for a, c in
                     zip(H[i, :j][::-1],
                         range(self.opening_cost, self.opening_cost + len(H[i, :j]) * (self.extention_cost),
                               self.extention_cost))]
        else:
            raise NotImplementedError(f'Gap type {gap_type} not implemented')
        max_cost = np.amax(costs)
        tuple = (max_cost, costs.index(max_cost) + 1)
        return tuple

    def _compute_scoring_matrix(self, *, ref: str, query: str):
        """
        Creates scoring matrix using the Smith Waterman algorithm with linear gap penalty.
        Keeps track of which value was computed using which neighbour in the traceback matrix.
        """
        H = np.zeros((len(ref) + 1, len(query) + 1), np.int)
        traceback_matrix = np.zeros((len(ref) + 1, len(query) + 1), dtype=object)
        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
            # DIAG: match or deletion
            match = H[i - 1, j - 1] + (
                self.match_score if ref[i - 1] == query[j - 1] else + self.mismatch_cost)
            # UP: deletion
            delete = self._compute_gap_costs(H, (i, j), 'D')
            # LEFT: insertion
            insert = self._compute_gap_costs(H, (i, j), 'I')
            scores = [match, delete[0], insert[0]]
            maximum = max(match, delete[0], insert[0], 0)
            H[i, j] = maximum
            traceback_matrix[i, j] = ['DIAG', ('UP', delete[1]), ('LEFT', insert[1])][
                scores.index(maximum)] if maximum in scores else 0
        self.traceback_matrix = traceback_matrix
        return H

    def _traceback(self, *, ref: str, query: str, loc: tuple):
        """
        Traces back the steps done for the computation of the scoring matrix.
        """
        (i, j) = loc
        am_i_done_yet = False  # end if encounter 0
        alignment = Alignment()
        alignment.end_coord = (i, j)
        alignment.score = self.score
        while not am_i_done_yet:
            if self.scoring_matrix[i, j] == 0:
                am_i_done_yet = True
            elif (self.traceback_matrix[i, j] == 'DIAG'):
                if query[j - 1] == ref[i - 1]:
                    alignment.prepend_to_cigar_string('=')
                else:
                    alignment.prepend_to_cigar_string('X')
                i -= 1
                j -= 1
            elif self.traceback_matrix[i, j][0] == 'UP':
                alignment.prepend_to_cigar_string('D', count=self.traceback_matrix[i, j][1])
                i -= (self.traceback_matrix[i, j][1])
            elif self.traceback_matrix[i, j][0] == 'LEFT':
                alignment.prepend_to_cigar_string('I', count=self.traceback_matrix[i, j][1])
                j -= (self.traceback_matrix[i, j][1])
        alignment.start_pos = i - 1  # -1 because 0-based
        alignment.start_coord = (i, j)
        yield alignment

    def __call__(self, *, ref: str, query: str):
        """
        Implements the Smith-Waterman alignment
        with linear gap penalty (same scores for opening and extending a gap)
        new score = max(
            match bonus  + prev score (i-1, j-1)
            substitution + prev score (i-1, j-1)
            gap penalty  + prev score (i-1, j)
            gap penalty  + prev score (i, j-1)
            0
        )

        Yields one alignment with maximal score per traceback,
        i.e. one for each last matching pair
        (assuming negative scores for mutation/indel)
        """
        self.scoring_matrix = self._compute_scoring_matrix(ref=ref, query=query)
        self.score = np.max(self.scoring_matrix)
        maxima = np.where(self.scoring_matrix == self.score)
        for loc in zip(*maxima):
            yield from self._traceback(ref=ref, query=query, loc=loc)


if __name__ == '__main__':
    ref = 'ATGGCCTC'
    query = 'ACGGCTC'
    aligner = SmithWaterman()
    for alignment in aligner(query=query, ref=ref):
        print(alignment.cigar_string)
        x, y, z = alignment.visualize(ref=ref, query=query)
        print(x)
        print(y)
        print(z)
        print(alignment.matching_subsegments())
        print(alignment.start_pos)
