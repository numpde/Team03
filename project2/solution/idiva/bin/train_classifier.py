from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from idiva.clf.utils import get_train_test
from idiva.clf.df import get_clinvar_clf_data
from idiva.clf.phenomenet import get_phenomenet
import typer
from enum import Enum
from keras.models import Sequential

"""
Options to consider:
     - encoding of the base-strings. Options: integer encoding, one-hot encoding, 
     only give information of start pos + length, 
     mutational signatures? (https://www.nature.com/articles/nature12477#Sec2)
Examples:
with integer encoding, training data looks like this:
           pos ref alt  label
0       949422   2   0      0
1       949422   2   0      0
2       949438   0   2      0
3       949438   0   2      0
4       949523   1   3      1

with only giving information of start pos + length:
           pos  label  length_var
0       949422      0           1
1       949422      0           1
2       949438      0           1
3       949438      0           1
4       949523      1           1

     - classifier type: random_forest, logistic regression
     
"""

CLASSIFIERS = {
    'log_reg': LogisticRegression(random_state=0, verbose=True, max_iter=1000),
    'random_forest': RandomForestClassifier(max_depth=2, random_state=0),
    'phenomenet': get_phenomenet
}


class BaseStringEncoding(Enum):
    integer = 'integer'
    base_string_length = 'base_string_length'


class WhichClassifier(Enum):
    log_reg = 'log_reg'
    random_forest = 'random_forest'
    phenomenet = 'phenomenet'


def main(base_string_encoding: BaseStringEncoding = typer.Argument(default='base_string_length',
                                                                   help='which encoding of the base-strings to use'),
         which_clf: WhichClassifier = typer.Argument(default='phenomenet', help='which classifier to use')):
    data_dir = Path(__file__).parent.parent.parent / 'data'

    # get the encoded data from clinvar
    clinvar_clf_data = get_clinvar_clf_data(data_dir, save_df=True, base_string_encoding=base_string_encoding.value)
    # split into train and validation sets
    train_data, train_labels, eval_data, eval_labels = get_train_test(clinvar_clf_data)

    if which_clf.value == 'phenomenet':
        clf: Sequential = CLASSIFIERS[which_clf.value](input_dim=train_data.shape[1]).fit(train_data, train_labels,
                                                                                          validation_data=(
                                                                                              eval_data, eval_labels),
                                                                                          batch_size=2500, verbose=2,
                                                                                          epochs=100)
    else:
        clf = CLASSIFIERS[which_clf.value].fit(train_data, train_labels)
        score = clf.score(eval_data, eval_labels)
        print(score)


if __name__ == '__main__':
    typer.run(main)
