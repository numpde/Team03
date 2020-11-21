# RA, 2020-11-20

import pandas
import scipy.stats


def v0_fisher(df: pandas.DataFrame, alternative='two-sided', oddsratio=False):
    """
    If oddsratio=True, compute the odds ratio.
    Otherwise compute p-value.

    From scipy doc:
        The calculated odds ratio is different from the one R uses.
        This scipy implementation returns the (more common)
        "unconditional Maximum Likelihood Estimate",
        while R uses the "conditional Maximum Likelihood Estimate".
        [https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html]
    """

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

        case0 = df[[ALTn_case]].to_numpy().sum(axis=1)
        case1 = df[set(case) - {ALTn_case}].to_numpy().sum(axis=1)

        ctrl0 = df[[ALTn_ctrl]].to_numpy().sum(axis=1)
        ctrl1 = df[set(ctrl) - {ALTn_ctrl}].to_numpy().sum(axis=1)

        matrixwise = (lambda arr: np.reshape(np.array(arr).T[:, :, np.newaxis], [-1, 2, 2]))

        assert np.array_equal(
            matrixwise([[10, 20, 30], [11, 21, 31], [12, 22, 32], [13, 23, 33]]),
            np.array([[[10, 11], [12, 13]], [[20, 21], [22, 23]], [[30, 31], [32, 33]]])
        )

        tests[F"ALT{n}_vs_Other"] = pd.Series(
            data=[
                # (log-odds ratio, p-value)
                (
                    (scipy.stats.fisher_exact(RxC, alternative=alternative))[1 - bool(oddsratio)]
                    if not np.isnan(RxC).any().any()
                    else
                    np.nan
                )
                for RxC in matrixwise([case0, case1, ctrl0, ctrl1])
            ],
            index=df.index,
        )

    return pd.DataFrame(data=tests)
