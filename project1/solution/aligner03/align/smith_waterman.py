import itertools
import numpy as np


class Smith_waterman_aligner(object):
    """
    Implementation of the Smith Waterman aligner algorithm
    """

    # todo:
    #  - integrate into workflow
    #  - implement better cost functions
    # https://elifesciences.org/articles/24284#fig1s1

    def __init__(self, ref, query, match_score, gap_cost, mismatch_cost):
        self.ref = ref
        self.query = query
        self.match_score = match_score
        self.gap_cost = gap_cost
        self.insertion_cost = gap_cost
        self.mismatch_cost = mismatch_cost
        self.cigar = []
        self.query_ = ''

    def create_scoring_matrix(self):
        """
        Creates scoring matrix
        """
        H = np.zeros((len(self.ref) + 1, len(self.query) + 1), np.int)

        for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
            match = H[i - 1, j - 1] + (
                self.match_score if self.ref[i - 1] == self.query[j - 1] else - self.mismatch_cost)
            delete = H[i - 1, j] - self.gap_cost
            insert = H[i, j - 1] - self.insertion_cost
            H[i, j] = max(match, delete, insert, 0)
        self.scoring_matrix = H
        return H

    def traceback(self, H, query_='', old_i=0, old_j=0):
        """
        Backtracing to find the optimal alignment
        returns: query_, start_pos
        """
        # flip H to get index of **last** occurrence of H.max() with np.argmax()
        H_flip = np.flip(np.flip(H, 0), 1)
        # find index of maximum in H_flip
        maxima = np.where(H_flip == H_flip.max())
        i_, j_ = (maxima[0].max(), maxima[1].max())
        # i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
        i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
        # distance makes only sense after first iteration when old values have been initialized
        dist_i = old_i - i if old_i - i > 0 else 1
        dist_j = old_j - j if old_j - j > 0 else 1
        # need to check if not traversing a 0 for ending condition
        for a in range(dist_i - 1):
            # if traversing 0 to find maxima, need to find maxima before 0
            if H_flip[a, a] == 0:
                i_ = j_ = a - 1
                i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
                dist_i = old_i - i
                dist_j = old_j - j
                query_ = self.query[j - 1:j - 1 + dist_j] + query_
                for idx_j, idx_i in zip(range(j - 1 + dist_j, j - 1, -1), range(i - 1 + dist_i, i - 1, -1)):
                    if self.ref[idx_i - 1] == self.query[idx_j - 1]:
                        # match
                        self.cigar.insert(0, '=')
                    else:
                        # mutation
                        self.cigar.insert(0, 'X')
                self.query_ = query_
                # zero is in the next cell, i.e. i-1
                self.start_pos = i - 1
                self.end_pos = i + len(query_) - 1
                return query_, self.start_pos
        if H[i, j] == 0:
            self.query_ = query_
            self.start_pos = i
            self.end_pos = i + 1 + len(query_)
            return query_, i

        if dist_i == dist_j:
            query_ = self.query[j - 1:j - 1 + dist_j] + query_
            for idx_j, idx_i in zip(range(j - 1 + dist_j, j - 1, -1), range(i - 1 + dist_i, i - 1, -1)):
                if self.ref[idx_i - 1] == self.query[idx_j - 1]:
                    # match
                    self.cigar.insert(0, '=')
                else:
                    # mutation
                    self.cigar.insert(0, 'X')
        elif dist_i > 1 and query_:
            # dist_i deletions
            query_ = self.query[j - 1] + '-' * (dist_i - 1) + query_
            for _ in range(dist_i - 1):
                self.cigar.insert(0, 'D')
            self.cigar.insert(0, '=')
        elif dist_j > 1 and query_ and not dist_i == dist_j:
            # dist_j insertions
            query_ = self.query[j - 1:j - 1 + dist_j] + query_
            for _ in range(dist_j - 1):
                self.cigar.insert(0, 'I')
            self.cigar.insert(0, '=')
        else:
            query_ = self.query[j - 1] + query_
            if self.ref[i - 1] == self.query[j - 1]:
                self.cigar.insert(0, '=')
            else:
                self.cigar.insert(0, 'X')

        assert len(self.cigar) == len(query_)
        return self.traceback(H[0:i, 0:j], query_, i, j)

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
                idx += 1

        return list

    def compress_cigar(self, cigar_list):
        counter = 1
        compressed = ''
        end = self.start_pos + len(cigar_list)
        for idx in range(len(cigar_list) - 1):
            if cigar_list[idx] == cigar_list[idx + 1]:
                counter += 1
            else:
                compressed += cigar_list[idx] + str(counter)
                counter = 1
        compressed += cigar_list[idx + 1] + str(counter)
        if self.start_pos:
            compressed = 'S' + str(self.start_pos) + compressed
        if (len(self.ref) - end):
            compressed += 'S' + str(len(self.ref) - end)
        return compressed

    def get_start_end(self):
        """
        Get start and end index of query in reference. 1-based!
        """
        # make reference and query uppercase
        self.ref, self.query = self.ref.upper(), self.query.upper()
        H = self.create_scoring_matrix()
        query_, pos = self.traceback(H)
        return pos + 1, pos + 1 + len(query_)

    def visualize(self, temp, read, cigar: str):
        """
        Taken from https://github.com/numpde/cbb/blob/master/20200610-CIGAR/cigar.py
        Example:
            (x, y, z) = visualize("CTAGCTACCGTCGTGCTA", "GCAACGACGATG", "S3=2X1=2D1=1X1=2I1=2S3")
        returns
            x = "TACCTAGCTACCGTCG-TGCTAGC"
            y = "       GCAAC-GACGATG    "
            z = "SSSSSS==X==D=X==I==SSSSS"
        """
        import re
        i = j = 0
        x = y = z = ""
        for (a, n) in re.findall(r"([=XIDS])([0-9]+)", cigar):
            n = int(n)
            z += a * n
            if (a in '=XDS'):
                x += temp[i:(i + n)]
                i += n
            if (a in '=XI'):
                y += read[j:(j + n)]
                j += n
            if (a == 'I'):
                x += "-" * n
            if (a == 'S'):
                y += " " * n
            if (a == 'D'):
                y += "-" * n
        return (x, y, z)


if __name__ == '__main__':
    test_waterman = Smith_waterman_aligner(ref='CTAGCTACCGTCGTGCTA', query='GCAACGACGATG',
                                           gap_cost=0, match_score=3,
                                           mismatch_cost=3)
    print(test_waterman.create_scoring_matrix())
    print(test_waterman.traceback(test_waterman.scoring_matrix))

    x, y, z = test_waterman.visualize(test_waterman.ref, test_waterman.query,
                                      test_waterman.compress_cigar(test_waterman.cigar))
    print(x)
    print(y)
    print(z)
