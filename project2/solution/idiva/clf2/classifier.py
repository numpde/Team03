import pandas as pd
from imblearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.svm import SVC

from idiva.dh.datahandler import DataHandler

class Classifier:
    """
    Object that classifies SNP variants into pathogenic or benign
    given a vcf file that contains CHR, POS, rsID, REF, ALT.
    The classifier is trained by default on clinvar.vcf.gz (GRCh37):
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
    """

    def __init__(self, clinvar_train: str = 'vcf_37'):
        self.model = None
        # Create model to train
        self.create_model()

        self.dataHandler = DataHandler()

        # create data for training
        self.x, self.labels = self.dataHandler.create_training_set(clinvar_train)

        # train model on training data
        self.train(self.x, self.labels)

    def create_model(self):
        """
        Returns the machine learning model
        """

        # TODO: Find a good model
        # TODO: What accuracy can we expect?

        selector = VarianceThreshold()
        scaler = StandardScaler()
        estimators = [
            ('mlp1', MLPClassifier(max_iter=1000)),
            ('rfc1', RandomForestClassifier())
        ]

        # classification
        classifier = StackingClassifier(estimators=estimators, final_estimator=SVC())

        steps = [('selector', selector),
                 ('scaler', scaler),
                 ('classification', classifier)]

        param_grid = {
            'classification__mlp1__alpha': np.logspace(start=-1, stop=1, num=5),

            'classification__rfc1__n_estimators': np.linspace(start=100, stop=500, num=5, dtype=np.int),

            'classification__final_estimator__C': np.logspace(start=-1, stop=1, num=5),
            'classification__final_estimator__kernel': ['poly']
        }

        pipeline = Pipeline(steps=steps)

        self.model = GridSearchCV(pipeline, param_grid, scoring='f1')

    def train(self, x_train: pd.DataFrame, labels: pd.DataFrame) -> None:
        """
        Fits the model to the given data
        """

        self.model.fit(x_train, labels.values.ravel())

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
