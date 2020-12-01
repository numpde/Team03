# RA, 2020-11-19

import scipy.stats
import pandas
import typing
import numpy


def chi2_test(df: pandas.DataFrame, cols: typing.Tuple[typing.List[str], typing.List[str]], add=0):
    """
    Chi-square test of independence of variables in a contingency table.
    The paired test is applied to each row of the dataframe `df`.
    One group of columns is cols[0] and the other is cols[1].

    Returns a series of p values.

    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2_contingency.html
    """

    NAs = df[cols[0] + cols[1]].isna().any(axis=1)

    return pandas.Series(
        data=(
            # chi2_contingency returns (chi2, p, dof, ex)
            (scipy.stats.chi2_contingency(RxC))[1] if not isNA else pandas.NA
            for (isNA, RxC) in zip(NAs, zip(df[cols[0]].to_numpy() + add, df[cols[1]].to_numpy() + add))
        ),
        index=df.index,
    )

