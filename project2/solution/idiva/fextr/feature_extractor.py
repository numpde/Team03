# LB DEZ 20
import os
import pickle
import warnings
from pathlib import Path
from typing import List

import numpy as np

import pandas as pd
import typing
from sklearn.dummy import DummyClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import Perceptron, SGDClassifier, PassiveAggressiveClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.utils import shuffle

from idiva import log
from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind

from tqdm import tqdm


class FeatureExtractor:
    """
    Given two vcf file this object trains a Perceptron classifier
    to be able to select the most important SNP's (GWAS with linear model)
    """

    def __init__(self, ctrl_vcf: str, case_vcf: str):
        self.clf, self.id = self.feature_extraction_chunks(ctrl_vcf, case_vcf)
        self.save_classifier()

    def save_classifier(self):
        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        filename = str(cache) + "/classifier.sav"

        pickle.dump(self.clf, open(filename, 'wb'))

    def get_extracted_variants(self) -> pd.DataFrame:
        """
        Returns the id's of the selected SNP's
        """
        selector = SelectFromModel(self.clf, prefit=True)

        return self.id[selector.get_support()]

    def get_reduced_dataframe(self) -> pd.DataFrame:
        """
        Returns reduced dataframe
        """
        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        from idiva.fextr import align

        with open(str(cache) + "/control_v2.vcf") as ctrl_vcf:
            ctrl_reader = ReadVCF(ctrl_vcf)
            with open(str(cache) + "/case_processed_v2.vcf") as case_vcf:
                case_reader = ReadVCF(case_vcf)
                dataframe = align(ctrl=ctrl_reader, case=case_reader)

        dataframe['ID'] = dataframe.ID_case.combine_first(dataframe.ID_ctrl)

        dataframe = dataframe[['CHROM', 'POS', 'ID', 'REF', 'ALT']]

        extracted = self.get_extracted_variants().values

        return dataframe.loc[extracted]

    @staticmethod
    def get_reduced_dataframe_from_saved_classifier() -> pd.DataFrame:
        """
        Returns reduced dataframe given that a classifier is stored as classifier.sav
        """
        clf = FeatureExtractor.get_saved_classifier()

        selector = SelectFromModel(clf, prefit=True)

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        from idiva.fextr import align

        with open(str(cache) + "/control_v2.vcf") as ctrl_vcf:
            ctrl_reader = ReadVCF(ctrl_vcf)
            with open(str(cache) + "/case_processed_v2.vcf") as case_vcf:
                case_reader = ReadVCF(case_vcf)
                dataframe = align(ctrl=ctrl_reader, case=case_reader)
                id = dataframe.index

        dataframe['ID'] = dataframe.ID_case.combine_first(dataframe.ID_ctrl)

        dataframe = dataframe[['CHROM', 'POS', 'ID', 'REF', 'ALT']]

        extracted = id[selector.get_support()].values

        return dataframe.loc[extracted]

    def feature_extraction_chunks(self, ctrl_vcf_file: str, case_vcf_file: str):
        """
        Returns a fitted Perceptron classifier for the given vcf files
        The classifier is trained in chunks where the chunks consist of a range of patient
        Therefore the classifier iterates columnwise over the vcf files
        The files are divided into equally many chunks and therefore the individual chunksize can differ
        """
        log.info("Fit linear classifier and reduce number of variants")

        clf = Perceptron()

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        # create unique index
        id = None

        with open(str(cache) + "/" + ctrl_vcf_file) as ctrl_vcf:
            ctrl_reader = ReadVCF(ctrl_vcf)
            with open(str(cache) + "/" + case_vcf_file) as case_vcf:
                case_reader = ReadVCF(case_vcf)
                dataframe = align(ctrl=ctrl_reader, case=case_reader)
                id = dataframe.index

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

                dataframe_ctrl['ID'] = dataframe_ctrl[['CHROM', 'POS', 'ALT']].apply(index_map, axis=1)

                dataframe_ctrl = dataframe_ctrl.drop(['CHROM', 'POS', 'ALT'], axis=1)

                dataframe_ctrl = dataframe_ctrl.set_index('ID')

                dataframe_ctrl = dataframe_ctrl.transpose()

                dataframe_ctrl = dataframe_ctrl.reindex(columns=id, fill_value=4)

                dataframe_case = pd.read_csv(fd_case, sep='\t', header=None,
                                             usecols=batch_columns_case,
                                             names=batch_names_case, converters=converter_dict_case)

                dataframe_case = dataframe_case.drop_duplicates(['CHROM', 'POS', 'ALT'], keep='first')

                dataframe_case['ID'] = dataframe_case[['CHROM', 'POS', 'ALT']].apply(index_map, axis=1)

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
        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        filename = str(cache) + "/classifier.sav"

        if os.path.exists(filename):
            loaded_model = pickle.load(open(filename, 'rb'))
        else:
            warnings.warn("no model saved")
            loaded_model = DummyClassifier()
        return loaded_model


def align(case: ReadVCF, ctrl: ReadVCF):
    """
    aligning case and control vcf file by joining on chrom, pos, ref and alt
    """
    from idiva.utils import seek_then_rewind

    dfs = {}
    for (k, vcf) in zip(['case', 'ctrl'], [case, ctrl]):
        with seek_then_rewind(vcf.fd, seek=vcf.dataline_start_pos) as fd:
            dfs[k] = pd.read_csv(fd, sep='\t', usecols=[0, 1, 2, 3, 4], header=None,
                                 names=["CHROM", "POS", "ID", "REF", "ALT"])
            dfs[k].index = dfs[k].index.rename(name="rowid")
            dfs[k] = dfs[k].reset_index().astype({'rowid': 'Int64'})

    dfs['case'] = dfs['case'].drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'], keep='first')
    dfs['ctrl'] = dfs['ctrl'].drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'], keep='first')

    df = join(case=dfs['case'], ctrl=dfs['ctrl'])

    df['CHROM'] = pd.to_numeric(df[['CHROM']].apply(translate_chrom, axis=1))

    df['CPA_ID'] = df[['CHROM', 'POS', 'ALT']].apply(index_map, axis=1)

    df = df.set_index('CPA_ID')

    # remove indels
    df = df[df['REF'].apply(lambda x: str(x) in ['A', 'C', 'G', 'T'])]
    df = df[df['ALT'].apply(lambda x: str(x) in ['A', 'C', 'G', 'T'])]

    return df


def join(case: pd.DataFrame, ctrl: pd.DataFrame) -> pd.DataFrame:
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


def index_map(chromposalt) -> int:
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


def translate_chrom(chrom: typing.Union[str, int]) -> int:
    """
    translate non integer chromosomes (X,Y & MT) to integers (23, 24 & 25)
    """

    if type(chrom) == pd.core.series.Series:
        chrom = chrom[0]

    if chrom == 'X':
        return 23
    elif chrom == 'Y':
        return 24
    elif chrom == 'MT':
        return 25
    else:
        return int(chrom)


if __name__ == '__main__':
    # fx = FeatureExtractor("control_v2.vcf", "case_processed_v2.vcf")
    # print("finished")

    clf = FeatureExtractor.get_saved_classifier()
    print(clf)
    print(type(clf))

    test = SelectFromModel(clf, prefit=True)
    print(test)
    print(test.get_support())
    print(sum(test.get_support()))

    print(FeatureExtractor.get_reduced_dataframe_from_saved_classifier())
    """

    import requests

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()

    file_path = os.path.join(cache, 'SIFT4G_Annotator2.jar')

    file_name = str(cache) + '/SIFT4G_Annotator2.jar'

    print(os.path.isfile(file_path))

    if not os.path.isfile(file_path):
        log.info("Downloading sift annotator (1.3MB)")
        results = requests.get('https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar')
        with open(file_name, 'wb') as f:
            f.write(results.content)

    print(os.path.getsize(file_path))

    dir_path = os.path.join(cache, 'GRCh37.74_new')

    dir_name_gz = str(cache) + '/GRCh37.74.zip'

    print(os.path.isdir(dir_path))

    if not os.path.isdir(dir_path):
        log.info("Downloading sift database (4.3 GB)")
        results = requests.get('https://sift.bii.a-star.edu.sg/sift4g/public//Homo_sapiens/GRCh37.74.zip')
        with open(dir_name_gz, 'wb') as f:
            f.write(results.content)

    import zipfile

    log.info("Unzip sift database (5.4 GB)")

    file = zipfile.ZipFile(dir_name_gz)
    file.extractall(path=dir_path)

    print(os.path.getsize(dir_path))

    # print(fx.get_extracted_variants())
    """