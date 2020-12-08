# HK, 2020-12-07

import idiva.io
from idiva import log
from idiva.clf.utils import get_trained_phenomenet
from idiva.dh.datahandler import DataHandler
from idiva.dh import datahandler
from idiva.clf2.classifier import Classifier


def phenomenet_classifier(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF) -> object:
    """
    Classifies the case-control df with a pretrained classifier.
    """
    from idiva.clf.df import c5_df

    log.info("Running the phenomenet classifier.")
    clf = Classifier()
    clf.model = get_trained_phenomenet()

    case_control = c5_df(case)
    results = clf.predict(case_control)
    case_control['class'] = results.values

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
