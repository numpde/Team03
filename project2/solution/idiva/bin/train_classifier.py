# HK, 2020-11-21
import os
from datetime import datetime

import pandas as pd
from idiva.clf.utils import TrainPhenomenetClinvardbSNPArgs
from idiva.clf2.classifier import Classifier
from dataclasses import dataclass
import typing


class Experiment:
    def __init__(self):
        if not os.path.exists('experiments_dataframe.csv'):
            self.experiments_dataframe = pd.DataFrame()
        else:
            self.experiments_dataframe = pd.read_csv('experiments_dataframe.csv')

        self.experiment_uid = self._create_experiment_uid()

    def _create_experiment_uid(self) -> str:
        dateTimeObj = datetime.now()
        return 'exp_' + dateTimeObj.strftime("%Y_%m_%d_%H_%M_%S_%f")

    def save_experiment(self, values: dict):
        values['experiment_uid'] = self.experiment_uid
        experiments_dataframe = self.experiments_dataframe.append(values, ignore_index=True)
        experiments_dataframe.to_csv('experiments_dataframe.csv', index=False)


def main():
    cv = Classifier()
    # # Create model to train
    # cv.create_model()
    # cv.x = cv.x
    # cv.labels = cv.labels
    # # train model on training data
    # cv.train(cv.x, cv.labels)
    # print(cv.model.cv_results_)
    exp = Experiment()
    # exp.save_experiment(cv.model.cv_results_)
    args = TrainPhenomenetClinvardbSNPArgs(weighted_loss=True, feature_list=['chrom', 'pos', 'var', 'label'])

    # train phenomenet
    epochs = 1
    batch_size = 2500
    history = cv.train_phenomenet_clinvar_dbSNP(epochs=epochs, batch_size=batch_size, args=args)
    values = {k: v[-1] for k, v in history.history.items()}
    values['model'] = 'phenomenet'
    values['epochs'] = epochs
    values['batch_size'] = batch_size
    exp.save_experiment(values)


if __name__ == '__main__':
    main()
