# RA, 2020-10-09


def unlist1(L):
    """
    Check that L has only one element at return it.
    """
    L = list(L)
    if not (len(L) == 1):
        raise ValueError(F"Expected an interable of length 1, got {len(L)}.")
    return L[0]
