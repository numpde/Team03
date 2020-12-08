# LB 23-11-2020


import typing

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.utils import class_weight

from idiva import log
from idiva.clf.phenomenet import Phenomenet, BestPhenomenet, get_model

from idiva.clf.utils import TrainPhenomenetArgs
from idiva.clf.utils import get_train_test
from idiva.dh.datahandler import DataHandler
from sklearn.model_selection import ParameterGrid
from idiva.clf.utils import get_trained_phenomenet


class Classifier:
    """
    Object that classifies SNP variants into pathogenic or benign
    given a vcf file that contains CHR, POS, rsID, REF, ALT.
    The classifier is trained by default on clinvar.vcf.gz (GRCh37):
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
    """

    def __init__(self, clinvar_train: str = 'vcf_37'):
        self.clinvar_train = clinvar_train
        self.model = None
        self.dataHandler = DataHandler()

        # create data for training
        self.x, self.labels = self.dataHandler.create_training_set()
        # self.phenomenet = Phenomenet(self.x.shape[1])

    def create_model(self, n_steps: int = 5, epochs: int = 100):
        """
        Returns the machine learning model
        """

        # TODO: Find a good model
        # TODO: What accuracy can we expect?
        # phenomenet = KerasClassifier(build_fn=self.phenomenet.get_phenomenet, batch_size=1, verbose=2, epochs=epochs)
        # phenomenet._estimator_type = "classifier"
        selector = VarianceThreshold()
        scaler = StandardScaler()
        estimators = [
            ('mlp1', MLPClassifier(max_iter=1000)),
            ('rfc1', RandomForestClassifier()),
            # ('phenomenet', phenomenet)
        ]

        # classification
        classifier = StackingClassifier(estimators=estimators, final_estimator=SVC())

        steps = [('selector', selector),
                 ('scaler', scaler),
                 ('classification', classifier)]

        param_grid = {
            'classification__mlp1__alpha': np.logspace(start=-1, stop=1, num=n_steps),

            'classification__rfc1__n_estimators': np.linspace(start=100, stop=500, num=n_steps, dtype=np.int),

            'classification__final_estimator__C': np.logspace(start=-1, stop=1, num=n_steps),
            'classification__final_estimator__kernel': ['poly']
        }

        pipeline = Pipeline(steps=steps)
        self.model = GridSearchCV(pipeline, param_grid, scoring='f1', verbose=2, n_jobs=-1)

    def get_phenomenet_data(self, args: TrainPhenomenetArgs):
        """
        HK, 2020-12-07
        """
        clf_data = self.dataHandler.get_phenomenet_training_data(args)

        if args.weighted_loss:
            weights = class_weight.compute_class_weight('balanced', np.unique(clf_data.label), clf_data.label)
            weights = {0: weights[0], 1: weights[1]}
        else:
            weights = None

        if args.feature_list:
            clf_data = clf_data[args.feature_list]

        weights: typing.Optional[typing.Iterable[float]]
        # split into train and validation sets
        train_data, train_labels, eval_data, eval_labels = get_train_test(
            clf_data, pipeline=Pipeline(steps=[('selector', VarianceThreshold()), ('scaler', StandardScaler())]))

        return train_data, train_labels, eval_data, eval_labels, weights

    def train_phenomenet(self, args: TrainPhenomenetArgs, phenomenet: typing.Optional[tf.keras.Sequential] = None) -> \
            typing.Tuple[
                tf.keras.callbacks.History, tf.keras.Sequential]:
        """
        Trains the phenomenet with the given parameters.
        HK, 2020-12-07
        """

        log.info(f'Training phenomenet with parameters: {args.__dict__}')
        train_data, train_labels, eval_data, eval_labels, weights = self.get_phenomenet_data(args)
        if args.which_phenomenet == 'basic':
            phenomenet = Phenomenet(train_data.shape[1])
            phenomenet = phenomenet.get_phenomenet()
        elif args.which_phenomenet == 'tuned':
            phenomenet = get_model(which=args.database)

        log.info(f'Training phenomenet for up to {args.epochs} epochs.')

        cb = tf.keras.callbacks.EarlyStopping(
            monitor='val_precision', min_delta=0, patience=args.early_stopping_patience, verbose=1, mode='max',
            baseline=None, restore_best_weights=True)

        history = phenomenet.fit(train_data, train_labels, validation_data=(
            eval_data, eval_labels),
                                 batch_size=args.batch_size, verbose=2,
                                 epochs=args.epochs, class_weight=weights, callbacks=[cb])

        return history, phenomenet

    def keras_tuner_rs(self, which_tuner: str, feature_list: typing.Optional[typing.Iterable[str]], database,
                       exp_name: str):
        """
        Launches keras tuner random search or hyperband optimization.
        HK, 2020-12-07
        """
        from idiva.clf.model_tuning.tuner import get_tuner

        args = TrainPhenomenetArgs(weighted_loss=True, database=database, feature_list=feature_list,
                                   early_stopping_patience=20)
        train_data, train_labels, eval_data, eval_labels, weights = self.get_phenomenet_data(args)

        cb = tf.keras.callbacks.EarlyStopping(
            monitor='val_precision', min_delta=0, patience=args.early_stopping_patience, verbose=1, mode='max',
            baseline=None, restore_best_weights=True)

        tuner = get_tuner(which_tuner, train_data.shape[1], exp_name)

        tuner.search_space_summary()

        tuner.search(x=train_data,
                     y=train_labels,
                     epochs=100,
                     validation_data=(eval_data, eval_labels), class_weight=weights,
                     batch_size=args.batch_size, callbacks=[cb], verbose=2)

        tuner.results_summary()
        best_model = tuner.get_best_models()[0]
        log.info('saving best tuner model.')
        best_model.save(f'best_tuner_model_{exp_name}')

    def train(self, x_train: pd.DataFrame, labels: pd.DataFrame) -> None:
        """
        Fits the model to the given data
        """

        self.model.fit(x_train.drop(['ID'], axis=1), labels.values.ravel())

    def predict(self, vcf_test1, vcf_test2=None) -> pd.DataFrame:
        """
        Returns predictions for the given vcf file
        """

        if vcf_test2 is not None:
            # fuse two files into one dataframe
            x_test = self.dataHandler.create_test_set_v2(case_vcf=vcf_test1, ctrl_vcf=vcf_test2)

        else:
            # create test features
            x_test = self.dataHandler.create_test_set(vcf_test1)

        # make predictions
        y_pred = self.model.predict(x_test)

        # create dataframe
        df = pd.DataFrame(y_pred)

        return df


def phenomenet_tuner_gird_search():
    search_space = {'params_clinvar_processd': {
        'feature_list': [None],
        'database': ['clinvar_processed'],
        'tuners': ['hyperband']
    },

        'params_clinvar_sbSNP': {
            'weighted_loss': [True],
            'feature_list': [['chrom', 'pos', 'var', 'label']],
            'database': ['clinvar_dbSNP'],
            'tuners': ['hyperband']
        }
    }

    cv = Classifier()
    for param_grid in ['params_clinvar_processd', 'params_clinvar_sbSNP']:
        for params in ParameterGrid(search_space[param_grid]):
            cv.keras_tuner_rs(which_tuner=params['tuners'], feature_list=params['feature_list'],
                              database=params['database'], exp_name='_'.join([params['tuners'], params['database']]))


if __name__ == '__main__':
    phenomenet_tuner_gird_search()
# # Create model to train
# cv.create_model(n_steps=1)
# cv.x = cv.x[:10]
# cv.labels = cv.labels[:10]
# # train model on training data
# cv.train(cv.x, cv.labels)
# print(cv.model.cv_results_)
