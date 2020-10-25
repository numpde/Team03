
from .naive import NaiveIndex

try:
    from .bw import BurrowsWheeler
    from .fm import FmIndex
except ModuleNotFoundError:
    import sys
    print(F"Warning: Import failed in {__file__}.", file=sys.stderr)
