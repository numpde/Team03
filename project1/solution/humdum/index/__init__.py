
from .naive import NaiveIndex

try:
    from .bw import BurrowsWheeler
    from .fm import FmIndex
except ModuleNotFoundError:
    pass
