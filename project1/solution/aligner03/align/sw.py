# HK, RA

import itertools
import numpy as np

from aligner03.utils import minidict
from aligner03.align import Alignment
from aligner03.align.alignment import prepend_to_cigar_string

default_mutation_costs = minidict({
    # Deletion
    'D': -2,
    # Insertion
    'I': -2,
    # Mutation
    'X': -1,
    # Match
    '=': 3,
})


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

    def _compute_scoring_matrix(self, *, ref: str, query: str):
        """
        Creates scoring matrix using the Smith Waterman algorithm with linear gap penalty.
        Keeps track of which value was computed using which neighbour in the traceback matrix.
        """
        H = np.zeros((len(ref) + 1, len(query) + 1), np.int)
        traceback_matrix = np.zeros((len(ref) + 1, len(query) + 1), np.int)
        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
            match = H[i - 1, j - 1] + (
                self.match_score if ref[i - 1] == query[j - 1] else + self.mismatch_cost)
            delete = H[i - 1, j] + self.gap_cost
            insert = H[i, j - 1] + self.insertion_cost
            scores = [match, delete, insert]
            maximum = max(match, delete, insert, 0)
            H[i, j] = maximum

            traceback_matrix[i, j] = scores.index(maximum) + 1 if maximum in scores else 0
        self.traceback_matrix = traceback_matrix
        return H

    def _traceback(self, *, ref: str, query: str, loc: tuple):
        """
        Traces back the steps done for the computation of the scoring matrix.
        """
        (i, j) = loc
        alignment = Alignment()
        alignment.end_coord = (i, j)
        alignment.score = self.score
        while 1:
            c = None

            if self.scoring_matrix[i, j] == 0:
                break
            if self.traceback_matrix[i, j] == 1:
                if query[j - 1] == ref[i - 1]:
                    c = '='
                else:
                    c = 'X'
                i -= 1
                j -= 1
            elif self.traceback_matrix[i, j] == 3:
                j -= 1
                c = 'D'
            elif self.traceback_matrix[i, j] == 2:
                i -= 1
                c = 'I'

            # This is inefficient: construct the string then compress [RA]
            alignment.cigar_string = prepend_to_cigar_string(c, alignment.cigar_string)

        alignment.start_pos = j + 1
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
