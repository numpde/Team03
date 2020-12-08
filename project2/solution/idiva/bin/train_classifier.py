# HK, 2020-11-21
import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from idiva import log
from idiva.clf.utils import TrainPhenomenetArgs
from idiva.clf2.classifier import Classifier
from sklearn.model_selection import ParameterGrid


class Experiment:
    def __init__(self):
        self.experiment_df = 'experiments_dataframe_new.csv'
        if not os.path.exists(self.experiment_df):
            self.experiments_dataframe = pd.DataFrame()
        else:
            self.experiments_dataframe = pd.read_csv(self.experiment_df)

        self.experiment_uid = self._create_experiment_uid()

    def _create_experiment_uid(self) -> str:
        dateTimeObj = datetime.now()
        return 'exp_' + dateTimeObj.strftime("%Y_%m_%d_%H_%M_%S_%f")

    def save_experiment(self, values: dict):
        log.info(f'Writing to experiment dataframe: {values}')
        values['experiment_uid'] = self.experiment_uid
        experiments_dataframe = self.experiments_dataframe.append(values, ignore_index=True)
        experiments_dataframe.to_csv(self.experiment_df, index=False)


search_space = {

    'params_clinvar_processd': {
        'feature_list': [None],
        'database': ['clinvar_processed'],
        'which_phenomenet': ['tuned']

    },

    'params_clinvar_sbSNP': {
        'feature_list': [['chrom', 'pos', 'var', 'label']],
        'database': ['clinvar_dbSNP'],
        'which_phenomenet': ['tuned']
    }
}


def train_phenomenet(args: TrainPhenomenetArgs):
    cv = Classifier()

    exp = Experiment()
    checkpoint_dir = Path(__file__).parent.parent / 'clf/checkpoints'
    checkpoint_dir.mkdir(exist_ok=True, parents=True)
    checkpoint_path = checkpoint_dir / exp.experiment_uid

    # train phenomenet
    history, model = cv.train_phenomenet(args=args)
    end_epoch = len(history.history['loss'])
    values = {k: v[-1] for k, v in history.history.items()}
    values['model'] = 'phenomenet'
    values['epochs'] = end_epoch
    values = {**values, **args.__dict__}
    exp.save_experiment(values)
    log.info(f'saving model to {checkpoint_path}_{end_epoch}')
    model.save(checkpoint_path)


if __name__ == '__main__':
    for param_grid in ['params_clinvar_processd', 'params_clinvar_sbSNP']:
        # for param_grid in ['params_clinvar_sbSNP']:
        for params in ParameterGrid(search_space[param_grid]):
            args = TrainPhenomenetArgs(**params)
            train_phenomenet(args=args)
