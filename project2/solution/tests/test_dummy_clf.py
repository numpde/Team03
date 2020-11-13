from itertools import product
from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split

from idiva.io.vcf import ReadVCF

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control.vcf",
    'case': BASE / "case_processed.vcf",
}


class NucEncoder:
    def __init__(self):
        self.n2i = {}
        self.i2n = {}
        for idx, key in enumerate(['A', 'C', 'G', 'T']):
            self.n2i[key] = idx
            self.i2n[idx] = key


def create_df(which_vcf: str) -> pd.DataFrame:
    """
    Creates a dataframe from a vcf file, containing the information of the POS, REF, ALT and SAMPLEs columns.
    The label is 1 for variants from the case vcf and 0 for variants from the control vcf.
    The nucleobase information from REF and ALT is transformed into and index.
    For each variant, the information of the SAMPLE column is included by inserting the number of
    occurences of values into the corresponding column. e.g. if the SAMPLES column contains [0|0, 1|0, 0|0],
    a 2 will be added to the "0|0" column of the dataframe and a 1 to the "1|0" column.

    Example:
            pos  ref  alt  label    0|0  0|1  0|2  1|0  1|1  1|2  2|0  2|1  2|2
        0  52.0  1.0  0.0    0.0  500.0  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
        1  56.0  1.0  3.0    0.0  498.0  1.0  NaN  1.0  NaN  NaN  NaN  NaN  NaN
        2  78.0  2.0  1.0    0.0  499.0  NaN  NaN  1.0  NaN  NaN  NaN  NaN  NaN
        3  80.0  2.0  0.0    0.0  497.0  1.0  NaN  2.0  NaN  NaN  NaN  NaN  NaN
        4  92.0  2.0  3.0    0.0  500.0  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    """
    df = pd.DataFrame(
        columns=['pos', 'ref', 'alt', 'label', *[f'{x}|{y}' for x, y in product(range(3), range(3))]], dtype='int32')
    nuc_encoder = NucEncoder()
    with PATHS[which_vcf].open(mode='r') as fd:
        for elem in ReadVCF(fd):
            line = {
                'pos': elem.pos,
                'ref': nuc_encoder.n2i[elem.ref],
                'alt': nuc_encoder.n2i[elem.alt],
                'label': 1 * (which_vcf == 'case'),
            }
            for sample, count in zip(*np.unique(elem.samples, return_counts=True)):
                line[sample] = int(count)
            df = df.append(line, ignore_index=True)
    return df


def get_datasets():
    df_ctrl = create_df('ctrl')
    df_case = create_df('case')
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
