from pathlib import Path

import numpy as np

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel

from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind


class FeatureExtractor:

    def __init__(self, ctrl_vcf: str, case_vcf: str):

        self.clf, self.id = self.feature_extraction_chunked(ctrl_vcf, case_vcf)

    def get_extracted_variants(self) -> pd.DataFrame:

        selector = SelectFromModel(self.clf, prefit=True)

        print(selector.get_support())

        return self.id[selector.get_support()]

    def feature_extraction_chunked(self, ctrl_vcf_file: str, case_vcf_file: str):
        """
        Returns a fitted RandomForestClassifier to the given vcf files
        The classifier is trained in chunks where the chunks consist of a range of patient
        Therefore the classifier iterates columnwise over the vcf files
        """
        rfc = RandomForestClassifier(n_estimators=100, warm_start=True)

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

        rfc = self.feature_extraction_file(ctrl_vcf_file, 'ctrl', id, rfc)
        rfc = self.feature_extraction_file(case_vcf_file, 'case', id, rfc)

        return rfc, id

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

            for idx, batch in enumerate(batches[:-1]):

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

                    print(dataframe)

                    dataframe = dataframe.reindex(columns=id, fill_value=4)

                    print(dataframe)

                    labels = None
                    if group == 'ctrl':
                        labels = np.ones(dataframe.shape[0])
                    else:
                        labels = np.zeros(dataframe.shape[0])

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
