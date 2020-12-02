# RA, 2020-12-02

import numpy as np
from unittest import TestCase


class TestFisher(TestCase):
    def test_pipfisher_vs_scipy(self):
        from fisher import pvalue as fisher_pip

        (a1, a2) = ([1, 10, 20], [20, 10, 1])
        # (c1, c2) = (["ALT0_case", "ALT1_case", "ALT2_case"], ["ALT0_ctrl", "ALT1_ctrl", "ALT2_ctrl"])
        # df = pd.DataFrame(index=[1], data=[a1 + a2], columns=[c1 + c2])

        c0 = fisher_pip(*[*[a1[0], a1[1] + a1[2]], *[a2[0], a2[1] + a2[2]]]).two_tail
        c1 = fisher_pip(*[*[a1[1], a1[0] + a1[2]], *[a2[1], a2[0] + a2[2]]]).two_tail
        c2 = fisher_pip(*[*[a1[2], a1[0] + a1[1]], *[a2[2], a2[0] + a2[1]]]).two_tail

        from idiva.stat import fisher_scipy

        (__, r0) = fisher_scipy([[a1[0], a1[1] + a1[2]], [a2[0], a2[1] + a2[2]]])
        (__, r1) = fisher_scipy([[a1[1], a1[0] + a1[2]], [a2[1], a2[0] + a2[2]]])
        (__, r2) = fisher_scipy([[a1[2], a1[0] + a1[1]], [a2[2], a2[0] + a2[1]]])

        self.assertEqual(F"{c0:.3e}", F"{r0:.3e}")
        self.assertEqual(F"{c1:.3e}", F"{r1:.3e}")
        self.assertEqual(F"{c2:.3e}", F"{r2:.3e}")

        with self.assertRaises(ValueError):
            fisher_pip(np.nan, np.nan, np.nan, np.nan)
