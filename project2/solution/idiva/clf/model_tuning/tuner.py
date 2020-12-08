import keras
import kerastuner
import tensorflow as tf
from kerastuner import HyperModel, Hyperband, RandomSearch
from tensorflow.keras import layers
from tensorflow.python.keras.layers import Dropout, Dense


class HyperPhenomenet(HyperModel):

    def __init__(self, input_dim):
        self.input_dim = input_dim

    def build(self, hp):
        model = keras.Sequential()
        model.add(layers.Dense(units=hp.Int('units_' + str(0),
                                            min_value=32,
                                            max_value=512,
                                            step=32),
                               activation='relu', input_dim=self.input_dim))
        for i in range(hp.Int('num_layers', 1, 20)):
            model.add(layers.Dense(units=hp.Int('units_' + str(i),
                                                min_value=32,
                                                max_value=512,
                                                step=32),
                                   activation='relu'))
            model.add(Dropout(rate=hp.Float('dropout_' + str(i), 0, 0.5, step=0.1, default=0.5)))
        model.add(Dense(1, kernel_initializer='uniform',
                        activation='sigmoid'))
        model.compile(
            optimizer=keras.optimizers.Adam(
                hp.Choice('learning_rate_Adam',
                          values=[1e-2, 1e-3, 1e-4])),
            loss='binary_crossentropy',

            metrics=['accuracy', 'AUC', tf.keras.metrics.Precision(name='precision'),
                     tf.keras.metrics.Recall()])
        return model


def get_tuner(which_tuner, input_shape, exp_name: str):
    tuners = {
        'hyperband': Hyperband(
            HyperPhenomenet(input_shape),
            objective=kerastuner.Objective("val_precision", direction="max"),
            directory=exp_name,
            project_name=exp_name,
            max_epochs=40
        ),
        'random_search': RandomSearch(
            HyperPhenomenet(input_shape),
            objective=kerastuner.Objective("val_precision", direction="max"),
            directory=exp_name,
            project_name=exp_name,
            max_trials=100)
    }
    return tuners[which_tuner]