# RA, 2020-12-04

# Template for a classifier interface to main

import idiva.io
from idiva import log


def placeholder(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):
    from idiva.clf.df import v0_df

    log.info("Running the placeholder classifier.")

    result = v0_df(case).assign(Dummy1=1, Dummy2=2)

    # The result should contain the columns CHROM, POS, ID
    # in the same order as `case`, but can have fewer rows

    class response:
        info = {
            'Dummy1': {'Number': '1', 'Type': 'Float', 'Description': '"Just the number 1"'},
            'Dummy2': {'Number': '1', 'Type': 'Float', 'Description': '"Just the number 2"'},
        }

        # This results in INFO lines as follows:
        # ##INFO=<ID=Dummy1,Number=1,Type=Float,Description="Just the number 1">
        # ##INFO=<ID=Dummy2,Number=1,Type=Float,Description="Just the number 2">

        df = result

    return response
