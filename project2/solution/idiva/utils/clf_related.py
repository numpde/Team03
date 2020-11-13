from itertools import product

import numpy as np
import pandas as pd
from joblib import load

from idiva.io import ReadVCF
from sklearn.dummy import DummyClassifier


def create_df(vcf_path: str, create_label=False, control: bool = False) -> pd.DataFrame:
    """
    Creates a dataframe from a vcf file, containing the information of the POS, REF, ALT and SAMPLEs columns.
    The nucleobase information from REF and ALT is transformed into and index.
    For each variant, the information of the SAMPLE column is included by inserting the number of
    occurrences of values into the corresponding column. e.g. if the SAMPLES column contains [0|0, 1|0, 0|0],
    a 2 will be added to the "0|0" column of the dataframe and a 1 to the "1|0" column.
    create_label (bool): indicates whether a label column is to be added to the dataframe.
        The label is 1 for variants from the case vcf and 0 for variants from the control vcf.
    control (bool): indicates if the vcf is from the control cohort


    Example:
            pos  ref  alt  label    0|0  0|1  0|2  1|0  1|1  1|2  2|0  2|1  2|2
        0  52.0  1.0  0.0    0.0  500.0  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
        1  56.0  1.0  3.0    0.0  498.0  1.0  NaN  1.0  NaN  NaN  NaN  NaN  NaN
        2  78.0  2.0  1.0    0.0  499.0  NaN  NaN  1.0  NaN  NaN  NaN  NaN  NaN
        3  80.0  2.0  0.0    0.0  497.0  1.0  NaN  2.0  NaN  NaN  NaN  NaN  NaN
        4  92.0  2.0  3.0    0.0  500.0  NaN  NaN  NaN  NaN  NaN  NaN  NaN  NaN
    """
    df = pd.DataFrame(
        columns=['pos', 'ref', 'alt', 'label', *[f'{x}|{y}' for x, y in product(range(3), range(3))]], dtype='int32')
    nuc_encoder = NucEncoder()
    with open(vcf_path, mode='r') as fd:
        for elem in ReadVCF(fd):
            line = {
                'pos': elem.pos,
                'ref': nuc_encoder.n2i[elem.ref],
                'alt': nuc_encoder.n2i[elem.alt],
            }
            if create_label:
                line['label'] = 1 * (not control)
            for sample, count in zip(*np.unique(elem.samples, return_counts=True)):
                line[sample] = int(count)
            df = df.append(line, ignore_index=True)
    return df


def get_clf(args):
    if args.which_clf == 'dummy':
        clf = DummyClassifier()
        x = np.ones((1, 12))
        y = np.ones((1, 1))
        clf.fit(x, y)
        # This is how you would load it otherwise:
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
        for idx, key in enumerate(['A', 'C', 'G', 'T', '-']):
            self.n2i[key] = idx
            self.i2n[idx] = key

    def encode(self, bases: str) -> int:
        """
        Encodes nucleobases into an integer
        """
        # todo this encoding makes no sense for longer bases
        if not bases:
            bases = '-'
        output = []
        for base in bases:
            output.append(str(self.n2i[base]))

        return int(''.join(output))
