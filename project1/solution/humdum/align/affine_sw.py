# HK, RA

import itertools
import numpy as np
from typing import List, Tuple
import re

from aligner03.align.sw import default_mutation_costs
from aligner03.align import Alignment

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
        Creates scoring matrix using the Smith Waterman algorithm with affine gap penalty.
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
        alignment._end_pair = (i, j)
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
        alignment._start_pos = i - 1  # -1 because 0-based
        alignment._start_pair = (i, j)
        yield alignment

    def __call__(self, *, ref: str, query: str):
        """
        Implements the Smith-Waterman alignment
        with affine gap penalty: implementation from Computational Biomedicine lecture 3, slide 26
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
        print(alignment.cigar)
        x, y, z = alignment.visualize(ref=ref, query=query)
        print(x)
        print(y)
        print(z)
        print(alignment.matching_subsegments())
        print(alignment._start_pos)
