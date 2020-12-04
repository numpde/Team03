from pathlib import Path

import numpy as np

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.utils import shuffle

from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind

from tqdm import tqdm


class FeatureExtractor:

    def __init__(self, ctrl_vcf: str, case_vcf: str):

        self.clf, self.id = self.feature_extraction_chunked(ctrl_vcf, case_vcf)

    def get_extracted_variants(self) -> pd.DataFrame:

        selector = SelectFromModel(self.clf, max_features=20000, prefit=True)

        print(selector.get_support())

        return self.id[selector.get_support()]

    def feature_extraction_chunked(self, ctrl_vcf_file: str, case_vcf_file: str):
        """
        Returns a fitted RandomForestClassifier to the given vcf files
        The classifier is trained in chunks where the chunks consist of a range of patient
        Therefore the classifier iterates columnwise over the vcf files
        """
        clf = RandomForestClassifier(n_estimators=100, warm_start=True)

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

        #clf = self.feature_extraction_file(ctrl_vcf_file, 'ctrl', id, clf)
        #clf = self.feature_extraction_file(case_vcf_file, 'case', id, clf)

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        def convert_strang(strang: str) -> int:
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

        def index_map(chromposalt) -> int:
            """
            uniquely map chrom, pos, alt to and index
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

                number_of_batches = int(max([np.ceil(len_ctrl/min_batch_size), np.ceil(len_case/min_batch_size)]))

                batch_size_ctrl = int(np.ceil(len_ctrl/number_of_batches))
                batch_size_case = int(np.ceil(len_case/number_of_batches))

                batches_ctrl = [i * batch_size_ctrl for i in range(number_of_batches)]
                batches_case = [i * batch_size_case for i in range(number_of_batches)]

                batches_ctrl.append(len_ctrl)
                batches_case.append(len_case)

                for idx in tqdm(range(number_of_batches), total=number_of_batches, postfix='feature selection'):

                    with seek_then_rewind(reader_ctrl.fd, seek=reader_ctrl.dataline_start_pos) as fd_ctrl:
                        with seek_then_rewind(reader_case.fd, seek=reader_case.dataline_start_pos) as fd_case:

                            batch_names_ctrl = names_ctrl[:3]
                            batch_names_case = names_case[:3]

                            batch_names_ctrl.extend(names_ctrl[batches_ctrl[idx]+3:batches_ctrl[idx + 1]+3])
                            batch_names_case.extend(names_case[batches_case[idx]+3:batches_case[idx + 1]+3])

                            batch_columns_ctrl = [0, 1, 4]
                            batch_columns_case = [0, 1, 4]

                            batch_columns_ctrl.extend(list(range(batches_ctrl[idx]+9, batches_ctrl[idx+1]+9)))
                            batch_columns_case.extend(list(range(batches_case[idx]+9, batches_case[idx+1]+9)))

                            converter_dict_ctrl = {}

                            for column in batch_names_ctrl:
                                if column not in ['CHROM', 'POS', 'ALT']:
                                    converter_dict_ctrl[column] = convert_strang

                            converter_dict_case = {}

                            for column in batch_names_case:
                                if column not in ['CHROM', 'POS', 'ALT']:
                                    converter_dict_case[column] = convert_strang

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

                            labels = np.ones(dataframe_ctrl.shape[0])
                            labels = np.append(labels, np.zeros(dataframe_case.shape[0]))

                            dataframe = dataframe_ctrl.append(dataframe_case)

                            dataframe, labels = shuffle(dataframe, labels, random_state=0)

                            clf.n_estimators += 100
                            clf.fit(dataframe, labels)

        return clf, id

    def feature_extraction_file(self, vcf_file: str, group: str, id: np.ndarray, clf):

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        def convert_strang(strang: str) -> int:
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

        def index_map(chromposalt) -> int:
            """
            uniquely map chrom, pos, alt to and index
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

        with open(str(cache) + "/" + vcf_file) as vcf:
            reader = ReadVCF(vcf)

            header = reader.header

            exclude = [2, 3, 5, 6, 7, 8]

            names = [i for idx, i in enumerate(header) if idx not in exclude]

            batch_size = 100

            number_batches = int(np.ceil((len(names) - 3) / batch_size))

            batches = [i * batch_size + 3 for i in range(number_batches)]
            batches.append(len(names))

            for idx, batch in tqdm(enumerate(batches[:-1]), total=len(batches)-1, postfix='train selection for '+group):

                with seek_then_rewind(reader.fd, seek=reader.dataline_start_pos) as fd:

                    batch_names = names[:3]
                    batch_names.extend(names[batches[idx]:batches[idx + 1]])

                    batch_columns = [0, 1, 4]
                    batch_columns.extend(list(range(batches[idx] + 6, batches[idx + 1] + 6)))

                    converter_dict = {}

                    for column in batch_names:
                        if column not in ['CHROM', 'POS', 'ALT']:
                            converter_dict[column] = convert_strang

                    dataframe = pd.read_csv(fd, sep='\t', header=None,
                                            usecols=batch_columns,
                                            names=batch_names, converters=converter_dict)

                    dataframe = dataframe.drop_duplicates(['CHROM', 'POS', 'ALT'], keep='first')

                    dataframe['ID'] = dataframe[['CHROM', 'POS', 'ALT']].apply(index_map, axis=1)

                    dataframe = dataframe.drop(['CHROM', 'POS', 'ALT'], axis=1)

                    dataframe = dataframe.set_index('ID')

                    dataframe = dataframe.transpose()

                    dataframe = dataframe.reindex(columns=id, fill_value=4)

                    labels = None
                    if group == 'ctrl':
                        labels = np.ones(dataframe.shape[0])
                    else:
                        labels = np.zeros(dataframe.shape[0])

                    print(clf)

                    clf.n_estimators += 100
                    clf.fit(dataframe, labels)

        return clf

    def align(self, case: ReadVCF, ctrl: ReadVCF):
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

        def index_map(chromposalt) -> int:
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

        df['ID'] = df[['CHROM', 'POS', 'ALT']].apply(index_map, axis=1)

        df = df.set_index('ID')

        return df

    def join(self, case: pd.DataFrame, ctrl: pd.DataFrame) -> pd.DataFrame:
        """
        Outer-join two dataframes on the columns CHROM, POS, ID.
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


if __name__ == '__main__':
    fx = FeatureExtractor("control_v2.vcf", "case_processed_v2.vcf")

    print(fx.get_extracted_variants())
