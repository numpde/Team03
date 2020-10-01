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
            self.query_ = query_
            return query_, i
        query_ = self.query[j - 1] + '-' + query_ if old_i - i > 1 else self.query[j - 1] + query_
        return self.traceback(H[0:i, 0:j], query_, i)

    def get_start_end(self):
        """
        Get start and end index of query in reference. 1-based!
        """
        # make reference and query uppercase
        self.ref, self.query = self.ref.upper(), self.query.upper()
        H = self.create_scoring_matrix()
        query_, pos = self.traceback(H)
        return pos + 1, pos + 1 + len(query_)

    def get_matching_blocks(self):
        """
        Get list of 1-based matching intervals
        Example:
            reference = ABCD
            query = ABD
            query_ = AB-D
            return: [[1, 2], [4, 4]]
        """
        list = []
        matching = False
        last_item = False
        self.ref, self.query = self.ref.upper(), self.query.upper()
        H = self.create_scoring_matrix()
        query_, pos = self.traceback(H)
        if query_ == self.ref[pos:pos + len(query_)]:
            # query_ matches exactly ref
            list.append([pos + 1, pos + len(query_)])
        else:
            idx = 1
            for a, b in zip(self.ref, query_):
                if idx == len(query_):
                    last_item = True
                if a == b:
                    if matching == False:
                        start_idx = idx
                        if last_item:
                            # single last query_ item corresponds to ref
                            list.append([pos + start_idx, pos + start_idx])
                    elif last_item:
                        list.append([pos + idx, pos + idx])
                    matching = True
                else:
                    if matching:
                        end_idx = idx - 1
                        list.append([pos + start_idx, pos + end_idx])
                    matching = False
                # print(idx, a, b, matching)
                idx += 1

        return list


if __name__ == '__main__':
    test_waterman = Smith_waterman_aligner(ref='ABCD', query='ABD', gap_cost=1, match_score=3, mismatch_cost=5)
    print(test_waterman.create_scoring_matrix())
    print(test_waterman.traceback(test_waterman.scoring_matrix))
    print(test_waterman.get_start_end())
    print(test_waterman.get_matching_blocks())
