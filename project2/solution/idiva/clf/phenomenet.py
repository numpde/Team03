# HK, 2020 - 11 - 21
import keras
import tensorflow as tf
from keras.layers import Dense
from keras.layers import Dropout
from keras.models import Sequential
from tensorflow import keras

""" 
taken from https://github.com/bio-ontology-research-group/phenomenet-vp/blob/master/dev/nn_final_training.py
"""


class Phenomenet:
    def __init__(self, input_dim: int):
        self.input_dim = input_dim

    def get_phenomenet(self):
        model = Sequential()
        model.add(Dense(67, input_dim=self.input_dim, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        model.add(Dense(32, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        model.add(Dense(256, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))

        adam = keras.optimizers.Adam(lr=0.001)
        loss = tf.keras.losses.BinaryCrossentropy()
        model.compile(loss=loss, optimizer=adam,
                      metrics=['accuracy', 'AUC', tf.keras.metrics.Precision(name='precision'),
                               tf.keras.metrics.Recall()])

        return model


def get_model(which: str):
    models = {
        'clinvar_dbSNP': '/home/hendrik/src/compbio/project2/solution/best_tuner_model_hyperband_clinvar_dbSNP',
        'clinvar_processed': '/home/hendrik/src/compbio/project2/solution/best_tuner_model_hyperband_clinvar_processed'
    }

    return tf.keras.models.load_model(models[which])


class BestPhenomenet:
    def __init__(self, input_dim: int):
        self.input_dim = input_dim

    def get_phenomenet(self):
        model = Sequential()
        # 0
        model.add(Dense(224, input_dim=self.input_dim, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        # 1
        model.add(Dense(256, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.1))
        # 2
        model.add(Dense(96, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        # 3
        model.add(Dense(192, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.3))
        # 4
        model.add(Dense(64, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        # 5
        model.add(Dense(352, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.4))
        # 6
        model.add(Dense(288, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.3))
        # 7
        model.add(Dense(448, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.3))
        # 8
        model.add(Dense(192, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        # 9
        model.add(Dense(448, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.2))
        # 10
        model.add(Dense(224, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0))
        # 11
        model.add(Dense(160, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.5))
        # 12
        model.add(Dense(512, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0))
        # 13
        model.add(Dense(256, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.1))
        # 14
        model.add(Dense(160, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.1))
        # 15
        model.add(Dense(288, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.3))
        # 16
        model.add(Dense(160, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.5))
        # 17
        model.add(Dense(32, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.5))
        # 18
        model.add(Dense(32, kernel_initializer='uniform', activation='relu'))
        model.add(Dropout(0.5))

        model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))

        adam = keras.optimizers.Adam(lr=0.0001)
        loss = tf.keras.losses.BinaryCrossentropy()
        model.compile(loss=loss, optimizer=adam,
                      metrics=['accuracy', 'AUC', tf.keras.metrics.Precision(name='precision'),
                               tf.keras.metrics.Recall()])

        return model
