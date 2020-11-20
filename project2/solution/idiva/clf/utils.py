from typing import Optional

import numpy as np
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import train_test_split


def get_clf(which_clf: str):
    if which_clf == 'dummy':
        clf = DummyClassifier()
        x = np.ones((1, 12))
        y = np.ones((1, 1))
        clf.fit(x, y)
        # This is how you would load it otherwise:
        # from joblib import load
        # clf = load('tests/data_for_tests/classifiers/dummy_clf.joblib')
    else:
        raise NotImplementedError
    return clf


class NucEncoder:
    """
    Contains two dictionaries. Nuc-to-index (n2i) contains the information to encode the nucleobase,
    and index-to-Nuc contains the information to decode the index into the nucleobase.
    A,C,G,T are the bases and "-" indicates an indel
    """

    def __init__(self):
        self.n2i = {}
        self.i2n = {}
        for idx, key in enumerate(['A', 'C', 'G', 'T', '-', 'N']):
            self.n2i[key] = idx
            self.i2n[idx] = key

    def encode(self, bases: Optional[str]) -> int:
        """
        Encodes nucleobases into an integer
        """
        # todo this encoding makes no sense for longer bases
        if not bases:
            bases = '-'
        output = [str(self.n2i[base]) for base in bases]
        return int(''.join(output))


def get_train_test(data):
    df_train, df_eval = train_test_split(data, test_size=0.2, shuffle=True)

    train_data = df_train.loc[:, df_train.columns != 'label'].to_numpy()
    train_labels = df_train['label'].to_numpy()
    eval_data = df_eval.loc[:, df_eval.columns != 'label'].to_numpy()
    eval_labels = df_eval['label'].to_numpy()

    return train_data, train_labels, eval_data, eval_labels
