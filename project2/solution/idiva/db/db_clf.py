# HK, 2020-12-05


import typing

import pandas as pd

import idiva.io
from idiva import log
from idiva.db import db


def classify(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF,
             case_control: typing.Optional[pd.DataFrame] = None) -> object:
    """
    Joins the case df and ctrl df on 'CHROM', 'POS', 'REF', 'ALT', 'ID' as a case-control df.
    Classifies the case-control df by querying the clinvar and dbSNP data.

    case_control: Possibility to pass the case-control dataframe directly
                  via the case_control input for testing purposes.
    """
    from idiva.clf.df import c5_df, join

    log.info("Running the database classifier.")
    if case_control is None:
        log.info("Joining case and control.")
        case_control = join(case=c5_df(case), ctrl=c5_df(ctrl), on=['CHROM', 'POS', 'REF', 'ALT', 'ID'])

    db_PosRefAlt = db.get_db_label_df()

    merge_on_PosRefAlt = case_control.merge(db_PosRefAlt, left_on=['POS', 'REF', 'ALT'], right_on=['pos', 'ref', 'alt'],
                                            how='left')
    merge_on_PosRefAlt['class'] = merge_on_PosRefAlt['class'].fillna(2)
    log.info(
        f"Found {len(merge_on_PosRefAlt) - merge_on_PosRefAlt.loc[merge_on_PosRefAlt['class'] == 2, 'class'].count()} "
        f"labels in databases")

    result = merge_on_PosRefAlt[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'class']]

    class response:
        id_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']

        info = {
            'class': {
                'Number': '.',
                'Type': 'Integer',
                'Description':
                    '"Number indicating to which class the variant belongs. '
                    '0 - Benign, 1 - Pathogenic, 2 - Unknown"'
            },
        }

        df = result

    assert set(response.id_cols).issubset(set(response.df.columns))
    assert set(response.info.keys()).issubset(set(response.df.columns))

    return response
