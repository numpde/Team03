import itertools
import numpy as np
from typing import List, Tuple
import re

mutation_costs = {
    # Deletion
    'D': -2,
    # Insertion
    'I': -2,
    # Mutation
    'X': -1,
    # Match
    '=': 3,
}


class Aligner:
    def __init__(self):
        self.cigar_string = ''
        self.start_pos = None
        self.start_coord = None  # (i,j) coordinates of start_pos
        self.end_coord = None  # (i,j) coordinates of end_pos
        self.score = None

    def __repr__(self):
        return F"Aligns at {self.start_pos} with CIGAR = {self.cigar_string} and score = {self.score}"

    def insert_to_cigar_string(self, symbol: str):
        if self.cigar_string and symbol == re.findall(r"[XIDS=]", self.cigar_string)[0]:
            current_count = re.findall(r"[0-9]+", self.cigar_string)[0]
            self.cigar_string = str(int(current_count) + 1) + self.cigar_string[len(current_count):]
        else:
            self.cigar_string = f'1{symbol}' + self.cigar_string

    def count_total(self) -> int:
        total = 0
        for count in re.findall(r"[0-9]+", self.cigar_string):
            total += int(count)
        return total

    def get_matching_blocks(self) -> List[Tuple]:
        cigar = self.cigar_string
        idx = self.start_pos
        matching_blocks = []
        for (n, a) in re.findall(r"([0-9]+)([XIDS=])", cigar):
            n = int(n)
            if a == '=' and idx == self.start_pos:
                matching_blocks.append((idx, idx + n))
            elif a == '=':
                matching_blocks.append((idx + 1, idx + n))
            idx += n

        return matching_blocks

    def visualize(self, ref, query):
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


class Smith_Waterman:
    def __init__(self, mutation_costs=mutation_costs):
        self.match_score = mutation_costs['=']
        self.mismatch_cost = mutation_costs['X']
        self.gap_cost = mutation_costs['D']
        self.insertion_cost = mutation_costs['I']

    def create_scoring_matrix(self, query, ref):
        """
        Creates scoring matrix
        """
        H = np.zeros((len(ref) + 1, len(query) + 1), np.int)

        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
            match = H[i - 1, j - 1] + (
                self.match_score if ref[i - 1] == query[j - 1] else + self.mismatch_cost)
            delete = H[i - 1, j] + self.gap_cost
            insert = H[i, j - 1] + self.insertion_cost
            H[i, j] = max(match, delete, insert, 0)
        return H

    def traceback(self, scoring_matrix, ref: str, query: str, max_pos: tuple):
        i, j = max_pos
        end_condition = False  # end if encounter 0
        aligner = Aligner()
        aligner.end_coord = (i, j)
        aligner.score = self.score
        while end_condition == False:
            diag_score = scoring_matrix[i - 1, j - 1]
            left_score = scoring_matrix[i - 1, j]
            upper_score = scoring_matrix[i, j - 1]
            best_neighbouring_score = max(diag_score, left_score, upper_score)
            if best_neighbouring_score == 0:
                end_condition = True
            if best_neighbouring_score == diag_score:
                if ref[j - 1] == query[i - 1]:
                    aligner.insert_to_cigar_string('=')
                else:
                    aligner.insert_to_cigar_string('X')
                i -= 1
                j -= 1
            elif best_neighbouring_score == upper_score:
                j -= 1
                aligner.insert_to_cigar_string('D')
            elif best_neighbouring_score == left_score:
                i -= 1
                aligner.insert_to_cigar_string('I')
        aligner.start_pos = j
        aligner.start_coord = (i, j)
        return aligner

    def __call__(self, ref, query):
        scoring_matrix = self.create_scoring_matrix(ref, query)
        self.score = np.max(scoring_matrix)
        maxima = np.where(scoring_matrix == self.score)
        i_maxima = maxima[0]
        j_maxima = maxima[1]
        for i, j in zip(i_maxima, j_maxima):
            yield self.traceback(scoring_matrix, max_pos=(i, j), ref=ref, query=query)


if __name__ == 'main':
    ref = 'SADWERWDFSDFS'
    query = 'AASDAWERSD'
    aligner = Smith_Waterman()
    for alignment in aligner(query=query, ref=ref):
        print(alignment.cigar_string)
        x, y, z = alignment.visualize(ref, query)
        print(x)
        print(y)
        print(z)
        print(alignment.get_matching_blocks())
