# HK, pre- 2020-10-09
# RA, 2020-10-09

import itertools
import numpy as np
from typing import List, Tuple
import re

from aligner03.utils import minidict

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


class Alignment:
    def __init__(self):
        self.cigar_string = ''
        self.start_pos = None
        self.start_coord = None  # (i,j) coordinates of start_pos
        self.end_coord = None  # (i,j) coordinates of end_pos
        self.score = None

    def __repr__(self):
        return F"Aligns at {self.start_pos} with CIGAR = {self.cigar_string} and score = {self.score}"

    def prepend_to_cigar_string(self, symbol: str):
        if self.cigar_string and symbol == re.findall(r"[XIDS=]", self.cigar_string)[0]:
            current_count = re.findall(r"[0-9]+", self.cigar_string)[0]
            self.cigar_string = str(int(current_count) + 1) + self.cigar_string[len(current_count):]
        else:
            self.cigar_string = f'1{symbol}' + self.cigar_string

    def count_total(self) -> int:
        """
        TODO: What does this count?
        """
        total = 0
        for count in re.findall(r"[0-9]+", self.cigar_string):
            total += int(count)
        return total

    def matching_subsegments(self) -> List[Tuple]:
        """
        TODO: Description
        0 or 1 based?
        """
        cigar = self.cigar_string
        idx = self.start_pos - 1
        segment = []
        for (n, a) in re.findall(r"([0-9]+)([XIDS=])", cigar):
            n = int(n)
            if a == '=':
                segment.append((idx + 1, idx + n))
            idx += n

        return segment

    def visualize(self, *, ref: str, query: str):
        """
        TODO: Description
        """
        cigar = self.cigar_string
        total = self.count_total()
        end_pos = self.start_pos + total

        query = query[self.start_coord[0]:self.end_coord[0]]
        if self.start_pos:
            cigar = str(self.start_pos) + 'S' + cigar
        if (len(ref) - end_pos):
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
    TODO: Usage example
    """
    def __init__(self, mutation_costs=default_mutation_costs):
        self.match_score = mutation_costs['=']
        self.mismatch_cost = mutation_costs['X']
        self.gap_cost = mutation_costs['D']
        self.insertion_cost = mutation_costs['I']

    def _compute_scoring_matrix(self, *, ref: str, query: str):
        """
        Creates scoring matrix
        TODO: Meaningful description
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
        TODO: Description
        """
        (i, j) = loc
        am_i_done_yet = False  # end if encounter 0
        alignment = Alignment()
        alignment.end_coord = (i, j)
        alignment.score = self.score
        while not am_i_done_yet:
            if self.scoring_matrix[i, j] == 0:
                am_i_done_yet = True
            if self.traceback_matrix[i, j] == 1:
                if query[j - 1] == ref[i - 1]:
                    alignment.prepend_to_cigar_string('=')
                else:
                    alignment.prepend_to_cigar_string('X')
                i -= 1
                j -= 1
            elif self.traceback_matrix[i, j] == 3:
                j -= 1
                alignment.prepend_to_cigar_string('D')
            elif self.traceback_matrix[i, j] == 2:
                i -= 1
                alignment.prepend_to_cigar_string('I')
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
