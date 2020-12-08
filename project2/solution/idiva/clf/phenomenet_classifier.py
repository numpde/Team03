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
    from idiva.fextr import align

    log.info("Running the phenomenet classifier.")
    clf = Classifier()
    clf.model = get_trained_phenomenet()

    results = clf.predict(case, ctrl)
    case_control = c5_df(case)
    log.info('Aligning case and control.')
    aligned = align(case, ctrl)
    aligned['phenomenet_class'] = results

    merge_on_PosRefAlt = case_control.merge(aligned, on=['POS', 'REF', 'ALT'], how='left')
    # filling all missing values with 2, for "unknown"
    merge_on_PosRefAlt['phenomenet_class'] = merge_on_PosRefAlt['phenomenet_class'].fillna(2)
    log.info(
        f"Classified {len(merge_on_PosRefAlt) - merge_on_PosRefAlt.loc[merge_on_PosRefAlt['phenomenet_class'] == 2, 'phenomenet_class'].count()} "
        f"labels.")

    class response:
        id_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']

        info = {
            'phenomenet_class': {'Number': '.',
                                 'Type': 'Integer',
                                 'Description': '"Number indicating to which class the variant belongs. '
                                                '0 - Benign, 1 - Pathogenic, 2 - Unknown"'
                                 },
        }

        df = case_control[[*id_cols, 'phenomenet_class']]

    assert set(response.id_cols).issubset(set(response.df.columns))
    assert set(response.info.keys()).issubset(set(response.df.columns))

    return response
