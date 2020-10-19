# HK, RA

import re
from typing import List, Tuple


class Alignment:
    def __init__(self):
        # CIGAR string
        self.cigar = ""
        # Alignment score
        self.score = None
        # Position of first match in the query:
        self._start_pos = None
        # (i,j) coordinates of start_pos:
        self._start_pair = None
        # (i,j) coordinates of end_pos
        self._end_pair = None

    @property
    def loc_in_query(self):
        """
        0-based.
        """
        return self._start_pair[1]

    @property
    def loc_in_ref(self):
        """
        0-based.
        """
        return self._start_pair[0]

    def __repr__(self):
        return F"(Alignment at {self.loc_in_query} => {self.loc_in_ref} with CIGAR `{self.cigar}` and score {self.score})"

    def prepend_to_cigar_string(self, c):
        raise NotImplementedError("This is deprecated, use the function outside of the class.")

    def matched_length(self) -> int:
        """
        Returns the total count of matches, deletions, mutations and insertions,
        i.e. the length of the expanded cigar string,
        i.e. the length of the matched segment on the reference string (incl. gaps).
        Example: for 2=1I3= returns 2+1+3 = 6.
        """
        total = 0
        for count in re.findall(r"[0-9]+", self.cigar):
            total += int(count)
        return total

    def matching_subsegments(self) -> List[Tuple]:
        """
        Returns a list of all matching 1-based subsegments of the query.
        Example:
            AAAAC
            AAABC
            returns: [(1,3), (5,5)]
        """
        cigar = self.cigar
        idx = self._start_pos - 1
        segment = []
        for (n, a) in re.findall(r"([0-9]+)([XIDS=])", cigar):
            n = int(n)
            if a == '=':
                segment.append((idx + 1, idx + n))
            idx += n

        return segment

    def visualize(self, *, ref: str, query: str):
        """
        Returns (x, y, z) where
            x: the reference
            y: the query
            z: the expanded cigar string
        for visualization.

        From
            https://github.com/numpde/cbb/tree/master/20200610-CIGAR
        """
        cigar = self.cigar
        total = self.matched_length()
        end_pos = self._start_pos + total

        query = query[self._start_pair[0]:self._end_pair[0]]
        if self._start_pos:
            cigar = str(self._start_pos) + 'S' + cigar
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

    def make_alignment_semilocal(self, query: str, mismatch_cost: int) -> None:
        """
        Makes the local alignment semi-local by adding a mismatch penalty to all mismatching elements before and
        after the alignment
        """
        # check if mismatch cost needs to be added in the beginning of alignment
        if self.loc_in_query:
            self.score += -self.loc_in_query * mismatch_cost
            self.cigar = f'{self.loc_in_query}X' + self.cigar
            # todo this only works if start_pos[0]>=start_pos[1]
            self._start_pair = (self._start_pair[0] - self.loc_in_query, self._start_pair[1] - self.loc_in_query)
        # check if mismatch cost needs to be added at the end of the alignment
        if len(query) - self._end_pair[1] > 0:
            mismatch_tail = len(query) - self._end_pair[1]
            self.score += -mismatch_tail * mismatch_cost
            self.cigar += f'{mismatch_tail}X'
            self._end_pair = (self._end_pair[0] + mismatch_tail, self._end_pair[1] + mismatch_tail)

    def compress_cigar(self, cigar_list: List):
        """
        Compresses a cigar list into a cigar string.
        Example:
            [=,=,X,I] -> 2=1X1I
        """
        counter = 1
        compressed = ''
        for idx in range(len(cigar_list)):
            if not idx == len(cigar_list) - 1 and cigar_list[idx] == cigar_list[idx + 1]:
                counter += 1
            else:
                compressed += str(counter) + cigar_list[idx]
                counter = 1
        self.cigar = compressed

