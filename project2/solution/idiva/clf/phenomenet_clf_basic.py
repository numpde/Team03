# HK, 2020-12-08

import idiva.io
from idiva import log
from idiva.clf.utils import get_trained_phenomenet
from idiva.dh import datahandler


def phenomenet_classifier_basic(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF) -> object:
    """
    Classifies the case-control df with a pretrained classifier.
    """
    from idiva.clf.df import c5_df

    log.info("Running the phenomenet classifier.")

    model = get_trained_phenomenet(which='exp_2020_12_08_11_05_12_615223')

    case_control = c5_df(case)

    case_control['var'] = case_control[['REF', 'ALT']].apply(lambda x: datahandler.MAPPING[x[0]][x[1]], axis=1)
    clf_data = case_control[['CHROM', 'POS', 'var']]

    # columns need to be in the same order than they were for the training of the classifier:
    clf_data = clf_data.reindex(sorted(clf_data.columns), axis=1)
    predictions = model.predict(clf_data.to_numpy().astype('float32'))
    case_control['phenom_class'] = predictions

    class response:
        id_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']

        info = {
            'phenom_class': {'Number': '.',
                             'Type': 'Integer',
                             'Description': '"Number predicted by the Phenomenet, indicating to which '
                                            'class the variant belongs. 0 - Benign, 1 - Pathogenic"'
                             },
        }

        df = case_control[[*id_cols, 'phenom_class']]

    assert set(response.id_cols).issubset(set(response.df.columns))
    assert set(response.info.keys()).issubset(set(response.df.columns))

    return response
