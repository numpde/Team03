# LB DEZ 20
import os
import pickle
import warnings
from pathlib import Path
from typing import List

import numpy as np

import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import Perceptron, SGDClassifier, PassiveAggressiveClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.utils import shuffle

from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind

from tqdm import tqdm


class FeatureExtractor:
    """
    Given two vcf file this object trains a Random Forest classifier to be able to select the most important SNP's
    """

    def __init__(self, ctrl_vcf: str, case_vcf: str):

        self.clf, self.id = self.feature_extraction_chunks(ctrl_vcf, case_vcf)

        filename = 'classifier.sav'
        pickle.dump(self.clf, open(filename, 'wb'))

    def get_extracted_variants(self) -> pd.DataFrame:
        """
        Returns the id's of the selected SNP's
        """
        selector = SelectFromModel(self.clf, max_features=20000, prefit=True)

        print(selector.get_support())

        return self.id[selector.get_support()]

    def feature_extraction_chunks(self, ctrl_vcf_file: str, case_vcf_file: str):
        """
        Returns a fitted RandomForestClassifier for the given vcf files
        The classifier is trained in chunks where the chunks consist of a range of patient
        Therefore the classifier iterates columnwise over the vcf files
        The files are divided into equally many chunks and therefore the individual chunksize can differ
        """
        # clf = RandomForestClassifier(n_estimators=10000, warm_start=True)
        clf = MultinomialNB()

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        # create unique index
        id = None

        with open(str(cache) + "/" + ctrl_vcf_file) as ctrl_vcf:
            ctrl_reader = ReadVCF(ctrl_vcf)
            with open(str(cache) + "/" + case_vcf_file) as case_vcf:
                case_reader = ReadVCF(case_vcf)
                dataframe = self.align(ctrl=ctrl_reader, case=case_reader)
                id = dataframe.index

        print(id)

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        with open(str(cache) + "/" + ctrl_vcf_file) as ctrl_vcf:
            with open(str(cache) + "/" + case_vcf_file) as case_vcf:

                reader_ctrl = ReadVCF(ctrl_vcf)
                reader_case = ReadVCF(case_vcf)

                header_ctrl = reader_ctrl.header

                header_case = reader_case.header

                exclude = [2, 3, 5, 6, 7, 8]

                names_ctrl = [i for idx, i in enumerate(header_ctrl) if idx not in exclude]
                names_case = [i for idx, i in enumerate(header_case) if idx not in exclude]

                len_ctrl = len(header_ctrl) - 9
                len_case = len(header_case) - 9

                min_batch_size = min([len_ctrl, len_case, 50])

                number_of_batches = int(max([np.ceil(len_ctrl / min_batch_size), np.ceil(len_case / min_batch_size)]))

                batch_size_ctrl = int(np.ceil(len_ctrl / number_of_batches))
                batch_size_case = int(np.ceil(len_case / number_of_batches))

                batches_ctrl = [i * batch_size_ctrl for i in range(number_of_batches)]
                batches_case = [i * batch_size_case for i in range(number_of_batches)]

                batches_ctrl.append(len_ctrl)
                batches_case.append(len_case)

                for idx in tqdm(range(number_of_batches), total=number_of_batches, postfix='feature selection'):

                    clf = self.feature_extraction_batch(reader_ctrl, reader_case, names_ctrl, names_case,
                                                        batches_ctrl, batches_case, idx, clf, id)

        return clf, id

    def feature_extraction_batch(self, reader_ctrl: ReadVCF, reader_case: ReadVCF, names_ctrl: List[str],
                                 names_case: List[str], batches_ctrl: List[int], batches_case: List[int], idx: int,
                                 clf, id: List[int]):
        """
        Returns a trained classifier on one batch
        loads from both files some patients (one batch) for training
        """
        with seek_then_rewind(reader_ctrl.fd, seek=reader_ctrl.dataline_start_pos) as fd_ctrl:
            with seek_then_rewind(reader_case.fd, seek=reader_case.dataline_start_pos) as fd_case:

                batch_names_ctrl = names_ctrl[:3]
                batch_names_case = names_case[:3]

                batch_names_ctrl.extend(names_ctrl[batches_ctrl[idx] + 3:batches_ctrl[idx + 1] + 3])
                batch_names_case.extend(names_case[batches_case[idx] + 3:batches_case[idx + 1] + 3])

                batch_columns_ctrl = [0, 1, 4]
                batch_columns_case = [0, 1, 4]

                batch_columns_ctrl.extend(list(range(batches_ctrl[idx] + 9, batches_ctrl[idx + 1] + 9)))
                batch_columns_case.extend(list(range(batches_case[idx] + 9, batches_case[idx + 1] + 9)))

                converter_dict_ctrl = {}

                for column in batch_names_ctrl:
                    if column not in ['CHROM', 'POS', 'ALT']:
                        converter_dict_ctrl[column] = self.convert_strang

                converter_dict_case = {}

                for column in batch_names_case:
                    if column not in ['CHROM', 'POS', 'ALT']:
                        converter_dict_case[column] = self.convert_strang

                dataframe_ctrl = pd.read_csv(fd_ctrl, sep='\t', header=None,
                                             usecols=batch_columns_ctrl,
                                             names=batch_names_ctrl, converters=converter_dict_ctrl)

                dataframe_ctrl = dataframe_ctrl.drop_duplicates(['CHROM', 'POS', 'ALT'], keep='first')

                dataframe_ctrl['ID'] = dataframe_ctrl[['CHROM', 'POS', 'ALT']].apply(self.index_map, axis=1)

                dataframe_ctrl = dataframe_ctrl.drop(['CHROM', 'POS', 'ALT'], axis=1)

                dataframe_ctrl = dataframe_ctrl.set_index('ID')

                dataframe_ctrl = dataframe_ctrl.transpose()

                dataframe_ctrl = dataframe_ctrl.reindex(columns=id, fill_value=4)

                dataframe_case = pd.read_csv(fd_case, sep='\t', header=None,
                                             usecols=batch_columns_case,
                                             names=batch_names_case, converters=converter_dict_case)

                dataframe_case = dataframe_case.drop_duplicates(['CHROM', 'POS', 'ALT'], keep='first')

                dataframe_case['ID'] = dataframe_case[['CHROM', 'POS', 'ALT']].apply(self.index_map, axis=1)

                dataframe_case = dataframe_case.drop(['CHROM', 'POS', 'ALT'], axis=1)

                dataframe_case = dataframe_case.set_index('ID')

                dataframe_case = dataframe_case.transpose()

                dataframe_case = dataframe_case.reindex(columns=id, fill_value=4)

                labels = np.zeros(dataframe_ctrl.shape[0])
                labels = np.append(labels, np.ones(dataframe_case.shape[0]))

                dataframe = dataframe_ctrl.append(dataframe_case)

                dataframe, labels = shuffle(dataframe, labels, random_state=0)

                """
                # for Random Forest Classifier
                clf.n_estimators += 10000
                clf.fit(dataframe, labels)
                """
                clf.partial_fit(dataframe, labels, classes=[0, 1])

        return clf

    def align(self, case: ReadVCF, ctrl: ReadVCF):
        """
        aligning case and control vcf file by joining on chrom, pos and alt
        """
        from idiva.utils import seek_then_rewind

        dfs = {}
        for (k, vcf) in zip(['case', 'ctrl'], [case, ctrl]):
            with seek_then_rewind(vcf.fd, seek=vcf.dataline_start_pos) as fd:
                dfs[k] = pd.read_csv(fd, sep='\t', usecols=[0, 1, 3, 4], header=None,
                                     names=["CHROM", "POS", "REF", "ALT"])
                dfs[k].index = dfs[k].index.rename(name="rowid")
                dfs[k] = dfs[k].reset_index().astype({'rowid': 'Int64'})

        dfs['case'] = dfs['case'].drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'], keep='first')

        dfs['ctrl'] = dfs['ctrl'].drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'], keep='first')

        df = self.join(case=dfs['case'], ctrl=dfs['ctrl'])

        df['ID'] = df[['CHROM', 'POS', 'ALT']].apply(self.index_map, axis=1)

        df = df.set_index('ID')

        return df

    def join(self, case: pd.DataFrame, ctrl: pd.DataFrame) -> pd.DataFrame:
        """
        Outer-join two dataframes on the columns CHROM, POS, ALT.
        Use the suffixes _case and _ctrl for the other ambiguous columns.

        RA, 2020-11-14
        LB, 2020-12-04 adapted
        """

        df = pd.merge_ordered(
            left=case, right=ctrl,
            suffixes=['_case', '_ctrl'],
            on=['CHROM', 'POS', 'REF', 'ALT'],
            how="outer",
        )

        return df

    def index_map(self, chromposalt) -> int:
        """
        Returns unique identifier by mapping chrom, pos & alt
        """
        chrom = chromposalt[0]
        pos = chromposalt[1]
        alt = chromposalt[2]

        if alt == "A":
            alt = 0
        elif alt == "C":
            alt = 1
        elif alt == "G":
            alt = 2
        elif alt == "T":
            alt = 3

        return chrom * 10000000000 + pos * 10 + alt

    def convert_strang(self, strang: str) -> int:
        """
        SNP as 0, 1, 2 for homozygous, heterozygous, and variant homozygous
        """
        if strang == "0|0":
            return 0
        elif strang == "0|1":
            return 1
        elif strang == "1|0":
            return 1
        elif strang == "1|1":
            return 2

        return np.nan

    @staticmethod
    def get_saved_classifier():
        """
        Returns the saved classifier if it exists
        otherwise a dummy classifier is returned
        """
        filename = 'classifier.sav'
        if os.path.exists(filename):
            loaded_model = pickle.load(open(filename, 'rb'))
        else:
            warnings.warn("no model saved")
            loaded_model = DummyClassifier()
        return loaded_model


if __name__ == '__main__':
    fx = FeatureExtractor("control_v2.vcf", "case_processed_v2.vcf")
    print("finished")
    clf = FeatureExtractor.get_saved_classifier()
    print(clf)
    print(type(clf))

    test = SelectFromModel(clf, prefit=True)
    print(test)
    print(test.get_support())
    print(sum(test.get_support()))

    # print(fx.get_extracted_variants())
