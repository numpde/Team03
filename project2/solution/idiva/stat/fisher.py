# RA, 2020-11-20

import pandas

# This is much faster than scipy.stat.fisher_exact
#
from fisher import pvalue as fisher_pip

def v0_fisher(df: pandas.DataFrame, alternative='two-sided', oddsratio=False):
    """
    This uses the pip package `fisher`.
    Only alternative='two-sided' is accepted.
    """

    from idiva import log

    if not (alternative == 'two-sided'): raise NotImplementedError
    if oddsratio: raise NotImplementedError

    import numpy as np
    import pandas as pd

    case = ["ALT0_case", "ALT1_case", "ALT2_case"]
    ctrl = ["ALT0_ctrl", "ALT1_ctrl", "ALT2_ctrl"]

    # Check existence of columns
    df[case]
    df[ctrl]

    tests = {}

    for n in [0, 1, 2]:
        # ALTn vs Rest
        (ALTn_case, ALTn_ctrl) = (F"ALT{n}_case", F"ALT{n}_ctrl")

        log.info(F"Computing p-values for ALT{n} vs Rest.")

        case0 = df[[ALTn_case]].to_numpy().sum(axis=1)
        case1 = df[set(case) - {ALTn_case}].to_numpy().sum(axis=1)

        ctrl0 = df[[ALTn_ctrl]].to_numpy().sum(axis=1)
        ctrl1 = df[set(ctrl) - {ALTn_ctrl}].to_numpy().sum(axis=1)

        matrixwise = (lambda arr: np.reshape(np.array(arr).T[:, :, np.newaxis], [-1, 2, 2]))

        assert np.array_equal(
            matrixwise([[10, 20, 30], [11, 21, 31], [12, 22, 32], [13, 23, 33]]),
            np.array([[[10, 11], [12, 13]], [[20, 21], [22, 23]], [[30, 31], [32, 33]]])
        )

        from tqdm import tqdm

        def fisher(RxC):
            if np.isnan(RxC).any().any():
                return np.nan
            else:
                return fisher_pip(*RxC.flatten()).two_tail

                # scipy alternative
                # return (scipy.stats.fisher_exact(RxC, alternative=alternative))[1 - bool(oddsratio)]

        from joblib import Parallel, delayed

        matrixwise = matrixwise([case0, case1, ctrl0, ctrl1])

        tests[F"ALT{n}_vs_Other"] = pd.Series(
            index=df.index,
            data=Parallel(n_jobs=7, prefer="processes")(
                delayed(fisher)(RxC)
                for RxC in tqdm(matrixwise)
            ),
        )

    return pd.DataFrame(data=tests)
