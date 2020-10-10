# RA, 2020-10-09


def unlist1(L):
    """
    Check that L has only one element at return it.
    """
    L = list(L)
    if not (len(L) == 1):
        raise ValueError(F"Expected an interable of length 1, got {len(L)}.")
    return L[0]


class minidict:
    """
    A slim version of a read-only dictionary.
    """

    def __init__(self, data: dict):
        self._data = dict(data)

    def __repr__(self):
        return repr(self._data)

    def __getitem__(self, item):
        return self._data[item]

    def __len__(self):
        return len(self._data)

    def keys(self):
        return self._data.keys()
