# LB 23-11-2020

import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.svm import SVC

from idiva.dh.datahandler import DataHandler
from idiva.clf.phenomenet import Phenomenet
from keras.wrappers.scikit_learn import KerasClassifier
from idiva.clf.utils import get_train_test
import keras


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
        self.x, self.labels = self.dataHandler.create_training_set(clinvar_train)
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

    def train_phenomenet(self, epochs=100, batch_size=2500) -> keras.callbacks.History:
        clinvar_clf_data = self.dataHandler.get_clinvar_clf_data(self.clinvar_train)
        # split into train and validation sets
        train_data, train_labels, eval_data, eval_labels = get_train_test(
            clinvar_clf_data[['CHROM', 'POS', 'VAR', 'label']],
            pipeline=Pipeline(steps=[('selector', VarianceThreshold()), ('scaler', StandardScaler())]))

        phenomenet = Phenomenet(train_data.shape[1])
        phenomenet = phenomenet.get_phenomenet()
        return phenomenet.fit(train_data, train_labels, validation_data=(
            eval_data, eval_labels),
                              batch_size=batch_size, verbose=2,
                              epochs=epochs)

    def train(self, x_train: pd.DataFrame, labels: pd.DataFrame) -> None:
        """
        Fits the model to the given data
        """

        self.model.fit(x_train.drop(['ID'], axis=1), labels.values.ravel())

    def predict(self, vcf_file_test: str, vcf_file_test2: str = None) -> pd.DataFrame:
        """
        Returns predictions for the given vcf file
        """

        if vcf_file_test2 is not None:
            # fuse two files into one dataframe
            x_test = self.dataHandler.create_test_set(vcf_file_test, vcf_file_test2)

        else:
            # create test features
            x_test = self.dataHandler.create_test_set(vcf_file_test)

        # make predictions
        y_pred = self.model.predict(x_test)

        # create dataframe
        df = pd.DataFrame(y_pred)

        return df


if __name__ == '__main__':
    cv = Classifier()
    # Create model to train
    cv.create_model(n_steps=1)
    cv.x = cv.x[:10]
    cv.labels = cv.labels[:10]
    # train model on training data
    cv.train(cv.x, cv.labels)
    print(cv.model.cv_results_)
    history = cv.train_phenomenet(epochs=3, batch_size=5)
    values = {k: v[-1] for k, v in history.history.items()}
    print(values)
