# RA, 2020-10-08

import typing


def find_all(pattern: str, template: str) -> typing.List[int]:
    """
    Find all occurrences of pattern string in template string.
    Returns a list of starting positions with 0- based indexing.
    """
    import re
    return [int(m.start()) for m in re.finditer("(?=" + pattern + ")", template)]


def forward(s):
    """
    Returns `s`.
    """
    return s


def reverse(s):
    """
    Genetic reverse-complement of `s`.
    """
    return ''.join({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}[c] for c in s[-1::-1])
