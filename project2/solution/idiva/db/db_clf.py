# HK, 2020-12-05

# Template for a classifier interface to main.
# Make a copy for your implementation.

import idiva.io
from idiva import log


def classify(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):
    from idiva.clf.df import c5_df, v0_df

    log.info("Running the placeholder 'classifier'.")

    # DUMMY result dataframe
    result = v0_df(case).assign(Dummy1=1, Dummy2=2.0)
    # result = c5_df(case).assign(Dummy1=1, Dummy2=2.0)

    # The result should contain the columns
    #   CHROM, POS, ID
    # or
    #   CHROM, POS, ALT
    # or
    #   CHROM, POS, REF, ALT
    # but can have fewer rows than `case`.
    #
    # Specify these columns in `id_cols`.

    class response:
        id_cols = ["CHROM", "POS", "ID"]

        info = {
            'Dummy1': {'Number': '1', 'Type': 'Integer', 'Description': '"Just the number 1"'},
            'Dummy2': {'Number': '1', 'Type': 'Float', 'Description': '"Just the number 2.0"'},
        }

        # This results in INFO lines as follows:
        # ##INFO=<ID=Dummy1,Number=1,Type=Integer,Description="Just the number 1">
        # ##INFO=<ID=Dummy2,Number=1,Type=Float,Description="Just the number 2.0">

        df = result

    assert set(response.id_cols).issubset(set(response.df.columns))
    assert set(response.info.keys()).issubset(set(response.df.columns))

    return response


def failure(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):
    pass
