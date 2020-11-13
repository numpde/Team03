import warnings
from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

from idiva.db import clinvar_open
from idiva.io import ReadVCF
from idiva.io.vcf import RawDataline
from idiva.utils import at_most_n
from idiva.clf.utils import create_df, NucEncoder

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control.vcf",
    'case': BASE / "case_processed.vcf",
}


def get_label_from_clinvar(clinvar_line: RawDataline):
    """
    Creates a label [0;1] from a clinvar vcf line. 1 indicates a sickness, 0 indicates no sickness.
    """
    temp_dict = {}
    for elem in clinvar_line.info.split(';'):
        k, v = elem.split('=')
        temp_dict[k] = v

    if 'CLNDN' not in temp_dict.keys():
        warnings.warn('CLNDN not specified', UserWarning)
        label = None
    else:
        if temp_dict['CLNDN'] == 'not_provided':
            label = 0
        else:
            label = 1
    return label


def create_df_clinvar():
    """
    Creates a dataframe from a clinvar vcf similar to "create_df".
    Example:
                pos ref alt  label
        0  865568.0   2   0    0.0
        1  865583.0   1   3    0.0
        2  865628.0   2   0    0.0
        3  865655.0   3   2    0.0
        4  865716.0   2   0    0.0
    """
    df = pd.DataFrame(
        columns=['pos', 'ref', 'alt', 'label'], dtype='int32')
    nuc_encoder = NucEncoder()
    with clinvar_open(which='vcf_37') as fd:
        vcf = ReadVCF(fd)
        for line in at_most_n(vcf.datalines, 5000):
            line = {
                'pos': line.pos,
                'ref': nuc_encoder.encode(line.ref),
                'alt': nuc_encoder.encode(line.alt),
                'label': get_label_from_clinvar(line),
            }
            # if "CLNDN" is not specified, label is none and current line is not added to the df
            if line['label']:
                try:
                    df = df.append(line, ignore_index=True)
                except Exception as e:
                    warnings.warn(f'{e}', UserWarning)
    return df


def get_datasets():
    df_ctrl = create_df(PATHS['ctrl'], create_label=True, control=True)
    df_case = create_df(PATHS['case'], create_label=True)
    df = df_ctrl.append(df_case, ignore_index=True).fillna(0).astype('int')
    df_train, df_eval = train_test_split(df, test_size=0.2, shuffle=True)
    train_data = df_train.loc[:, df_train.columns != 'label'].to_numpy()
    train_labels = df_train['label'].to_numpy()
    eval_data = df_eval.loc[:, df_eval.columns != 'label'].to_numpy()
    eval_labels = df_eval['label'].to_numpy()
    return train_data, train_labels, eval_data, eval_labels


class TestClf(TestCase):
    def test_dummy_classifier(self):
        dummy_clf = DummyClassifier(strategy="most_frequent")
        train_data, train_labels, eval_data, eval_labels = get_datasets()
        dummy_clf.fit(train_data, train_labels)
        dummy_clf.predict(train_data)
        score = dummy_clf.score(train_data, train_labels)
        self.assertEqual(0.5, np.round(score, 1))

    def test_logreg(self):
        train_data, train_labels, eval_data, eval_labels = get_datasets()
        clf = LogisticRegression(random_state=0).fit(train_data, train_labels)
        score = clf.score(eval_data, eval_labels)
        self.assertIsNotNone(score)

    def test_create_df(self):
        df_ctrl = create_df(PATHS['ctrl'], create_label=True, control=True)
        df_case = create_df(PATHS['case'], create_label=True)
        df = df_ctrl.append(df_case, ignore_index=True).fillna(0).astype('int')
        self.assertTrue(len(df))

    def test_create_df_clinvar(self):
        df = create_df_clinvar()
        self.assertTrue(len(df))
