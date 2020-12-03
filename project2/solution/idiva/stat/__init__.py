# RA, 2020-12-02

# RuntimeWarning: numpy.ufunc size changed
# https://github.com/numpy/numpy/issues/14920#issuecomment-554672523
import warnings
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

from .fisher import v0_fisher
from .chi2 import chi2_test
from scipy.stats import fisher_exact as fisher_scipy
