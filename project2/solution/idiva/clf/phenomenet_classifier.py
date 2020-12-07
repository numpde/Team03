# HK, 2020-12-07

import idiva.io
from idiva import log
from idiva.clf.utils import get_trained_phenomenet
from idiva.dh import datahandler


def phenomenet_classifier(*, case: idiva.io.ReadVCF) -> object:
    """
    Classifies the case-control df by querying the clinvar and dbSNP data.
    """
    from idiva.clf.df import c5_df

    log.info("Running the phenomenet classifier.")

    model = get_trained_phenomenet()

    case_control = c5_df(case)

    case_control['var'] = case_control[['REF', 'ALT']].apply(lambda x: datahandler.MAPPING[x[0]][x[1]], axis=1)
    clf_data = case_control[['CHROM', 'POS', 'var']]
    # ----------------------------------------------------------------------------
    # todo get CASS scores
    for col in ['CADD_PHRED', 'CADD_SUCC', 'SIFT_SCORE', 'SIFT_SUCC']:
        clf_data[col] = 1
    # ----------------------------------------------------------------------------
    # todo columns need to be in the right order

    predictions = model.predict(clf_data.to_numpy().astype('float32'))
    case_control['class'] = predictions

    class response:
        id_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']

        info = {
            'class': {'Number': '.',
                      'Type': 'Integer',
                      'Description': '"Number indicating to which class the variant belongs. '
                                     '0 - Benign, 1 - Pathogenic, 2 - Unknown"'
                      },
        }

        df = case_control[[*id_cols, 'class']]

    assert set(response.id_cols).issubset(set(response.df.columns))
    assert set(response.info.keys()).issubset(set(response.df.columns))

    return response
