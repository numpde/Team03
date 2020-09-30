import itertools
import numpy as np


class Smith_waterman_aligner(object):
    """
    Toy implementation from https://tiefenauer.github.io/blog/smith-waterman/ to play around
    """

    # todo:
    #  - integrate into workflow
    #  - implement better cost functions

    def __init__(self, ref, query, match_score, gap_cost, mismatch_cost):
        self.ref = ref
        self.query = query
        self.match_score = match_score
        self.gap_cost = gap_cost
        self.mismatch_cost = mismatch_cost

    def create_scoring_matrix(self):
        """
        Creates scoring matrix
        """
        H = np.zeros((len(self.ref) + 1, len(self.query) + 1), np.int)

        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
            match = H[i - 1, j - 1] + (
                self.match_score if self.ref[i - 1] == self.query[j - 1] else - self.mismatch_cost)
            delete = H[i - 1, j] - self.gap_cost
            insert = H[i, j - 1] - self.gap_cost
            H[i, j] = max(match, delete, insert, 0)
        self.scoring_matrix = H
        return H

    def traceback(self, H, query_='', old_i=0):
        """
        Backtracing to find the optimal alignment
        """
        # flip H to get index of **last** occurrence of H.max() with np.argmax()
        H_flip = np.flip(np.flip(H, 0), 1)
        i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
        i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
        if H[i, j] == 0:
            return query_, j
        query_ = self.query[j - 1] + '-' + query_ if old_i - i > 1 else self.query[j - 1] + query_
        return self.traceback(H[0:i, 0:j], query_, i)

    def get_start_end(self):
        # make reference and query uppercase
        self.ref, self.query = self.ref.upper(), self.query.upper()
        H = self.create_scoring_matrix()
        b_, pos = self.traceback(H)
        return pos, pos + len(b_)


if __name__ == '__main__':
    test_waterman = Smith_waterman_aligner(ref='ATGGCCTC', query='ACGGCTC', gap_cost=4, match_score=1, mismatch_cost=3)
    print(test_waterman.create_scoring_matrix())
    print(test_waterman.traceback(test_waterman.scoring_matrix))
    print(test_waterman.get_start_end())
