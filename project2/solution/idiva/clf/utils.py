from typing import Optional

import numpy as np
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from idiva import log
import typing
from dataclasses import dataclass
from pathlib import Path
import os
import tensorflow as tf


@dataclass
class TrainPhenomenetArgs:
    weighted_loss: bool
    feature_list: typing.Optional[typing.Iterable[str]]
    database: str
    epochs: int = 1000
    batch_size: int = 500
    early_stopping_patience: int = 100


def get_trained_phenomenet():
    """
    Downloads trained phenomenet from polybox if not found in download folder
    """
    base = Path(__file__).parent.parent / 'download'

    exp_str = 'exp_2020_12_07_14_53_54_611564'
    model_path = base / exp_str
    if not model_path.exists():
        command = f'wget -O {base / exp_str}.tar.gz https://polybox.ethz.ch/index.php/s/bSwajA85feMAFiq/download'
        log.info(f'Downloading phenomenet with command {command}')
        os.system(command)
        command = f'tar -zxvf {base / exp_str}.tar.gz -C {base}'
        log.info(f'Extracting model with command {command}')
        os.system(command)
        log.info('removing compressed folder')
        os.system(f'rm {base / exp_str}.tar.gz')
    log.info('Loading trained phenomenet.')
    return tf.keras.models.load_model(model_path)


def get_clf(which_clf: str):
    """
    HK, 2020-11-21
    """
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
    HK, 2020-11-12
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


def get_train_test(data, pipeline: Pipeline = None):
    """
    HK, 2020-11-12
    """
    log.info('Splitting and preprocessing data.')
    df_train, df_eval = train_test_split(data, test_size=0.2, shuffle=True)

    train_data = df_train.loc[:, df_train.columns != 'label'].to_numpy()
    train_labels = df_train['label'].to_numpy()
    eval_data = df_eval.loc[:, df_eval.columns != 'label'].to_numpy()
    eval_labels = df_eval['label'].to_numpy()
    pipeline.fit(train_data, train_labels)
    train_data = pipeline.transform(train_data)
    eval_data = pipeline.transform(eval_data)
    return train_data, train_labels, eval_data, eval_labels
