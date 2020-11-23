# LB 23-11-2020

import typing

import numpy as np
import pandas as pd
from tqdm import tqdm

from idiva.db.ncbi_scraper import get_functional_consequence_from_SNP
from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind


class DataHandler:
    URL = {
        'vcf_37': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
        'vcf_38': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    }

    """
    A -> C = 0
    A -> G = 1
    A -> T = 2

    C -> A = 3
    C -> G = 4
    C -> T = 5

    G -> A = 6
    G -> C = 7
    G -> T = 8

    T -> A = 9
    T -> C = 10
    T -> G = 11
    """

    mapping = {'A': {'C': 1, 'G': 2, 'T': 2},
               'C': {'A': 3, 'G': 4, 'T': 5},
               'G': {'A': 6, 'C': 7, 'T': 8},
               'T': {'A': 9, 'C': 10, 'G': 11}}

    INIT_COLS = ["CHROM", "POS", "ID", "REF", "ALT"]
    CLINVAR_COLS = INIT_COLS
    CLINVAR_COLS_IDX = [0, 1, 2, 3, 4]

    def get_clf_datalines(self, df_clinvar: pd.DataFrame):
        """
        HK, 2020-11-21
        """
        for idx, row in tqdm(df_clinvar.iterrows(), total=len(df_clinvar), postfix='iterating df_clinvar'):
            if (str(row.ref) in ['A', 'C', 'G', 'T']) and (str(row.alt) in ['A', 'C', 'G', 'T']):
                p_succes, p_score = self.get_polyphen2_score(row.id)
                s_succes, s_score = self.get_sift_score(row.id)
                c_succes, c_score = self.get_cadd_score(row.id)

                line = {
                    'ID': row.id,
                    'CHROM': row.chrom,
                    'POS': row['pos'],
                    'VAR': self.mapping[row.ref][row.alt],
                    'label': 1 if row['CLNSIG'] == 'Pathogenic' else 0,
                    'PS': p_score,
                    'PB': p_succes,
                    'SS': s_score,
                    'SB': s_succes,
                    'CS': c_score,
                    'CB': c_succes,
                }

                yield line

    def df_clinvar_to_clf_data(self, df_clinvar: pd.DataFrame) -> pd.DataFrame:
        """
        HK, 2020-11-21
        """
        return pd.DataFrame(data=self.get_clf_datalines(df_clinvar))

    def get_clinvar_clf_data(self, clinvar_file: str = 'vcf_37') -> pd.DataFrame:
        """
        Loads clinvar_clf_data suitable for a classifier.

        HK, 2020-11-22
        RA, 2020-11-22
        """

        from idiva.io import cache_df

        which = 'vcf_37'

        def maker_clinvar(which) -> pd.DataFrame:
            from idiva.db import clinvar_open
            from idiva.io import ReadVCF
            from idiva.db.clinvar import clinvar_to_df

            with clinvar_open(which=which) as fd:
                return clinvar_to_df(ReadVCF(fd))

        df_clinvar = cache_df(name=("clinvar_" + which), key=[], df_maker=maker_clinvar)
        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]

        return cache_df(name="clinvar_clf_data", key=["v01"], df_maker=self.df_clinvar_to_clf_data,
                        df=df_clinvar_reduced)

    def create_training_set(self, clinvar_file: str = 'vcf_37') -> typing.Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Returns training features and corresponding labels given a clinvar vcf file
        """

        # create training set containg
        # CHROM, POS, VAR, Polyphen2 score & success, sift score & success, cadd score & success
        clinvar_clf_data = self.get_clinvar_clf_data()

        x_train = clinvar_clf_data.loc[:, clinvar_clf_data.columns != 'label']
        y_train = clinvar_clf_data.loc[:, clinvar_clf_data.columns == 'label']
        y_train.set_index(x_train.index)

        return x_train, y_train

    def create_test_set(self, vcf_file: str, vcf_file2: str = None) -> pd.DataFrame:

        frame = self.translate_vcf(vcf_file)

        if vcf_file2 is not None:
            baseframe2 = self.translate_vcf(vcf_file2)

            # merge both frames into one
            frame = pd.concat([frame, baseframe2])
            frame.drop_duplicates()

        return frame

    def get_labels(self, x: pd.DataFrame) -> pd.DataFrame:
        """
        Returns a dataframe containing the true labels for a given vcf file
        """
        # get Id
        ids = x['Id']

        labels = np.zeros(ids.size)

        for id, rs in enumerate(ids):
            # TODO Interpretation function for functional consequence => benign or pathogenic ?
            labels[id] = get_functional_consequence_from_SNP(rs)

        return pd.DataFrame(data=labels, columns=['label'])

    def translate_vcf(self, vcf_file: str) -> pd.DataFrame:
        """
        Returns a dataframe that contains the following features from a vcf file
        CHROM, POS, ID, VAR
        """

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        with open(str(cache) + "/" + vcf_file) as vcf:
            reader = ReadVCF(vcf)

            with seek_then_rewind(reader.fd, seek=reader.dataline_start_pos) as fd:

                dataframe = pd.read_csv(fd, sep='\t', usecols=range(len(DataHandler.INIT_COLS)), header=None,
                                        names=DataHandler.INIT_COLS,
                                        dtype={'CHROM': np.int, 'POS': np.int, 'ID': np.str, 'REF': np.str,
                                               'ALT': np.str})

                # Check if ALT contains only one value or several values seperated by ','
                assert (len([uni for uni in dataframe['ALT'].unique().tolist() if ',' in uni]) == 0)

                # store only SNP variants
                dataframe = dataframe[dataframe['REF'].apply(lambda x: set([x]).issubset({'A', 'C', 'G', 'T'}))]
                dataframe = dataframe[dataframe['ALT'].apply(lambda x: set([x]).issubset({'A', 'C', 'G', 'T'}))]

                # Check if only SNP
                for ref in dataframe['REF']:
                    assert (len(ref) == 1)

                for alt in dataframe['ALT']:
                    assert (len(alt) == 1)

                assert (set(dataframe['REF'].unique().tolist()).issubset({'A', 'C', 'G', 'T'}))
                assert (set(dataframe['ALT'].unique().tolist()).issubset({'A', 'C', 'G', 'T'}))

        dataframe['CHROM'] = pd.to_numeric(dataframe[['CHROM']].apply(self.translate_chrom, axis=1))

        dataframe = self.encode_ref_alt(dataframe)

        dataframe.drop_duplicates()

        # TODO:        same CHROM POS and rsID but not same REF & ALT
        #              same CHROM rsID REF ALT but not same POS
        #              => rsIDs are not completely unique !
        #              Ignore rsID (Kjong Nov 23)
        """
        
        print(len(dataframe['ID'].unique().tolist()))
        print(len(dataframe['ID'].tolist()))

                 CHROM       POS           ID REF ALT  VAR
        56638       17   1649616  rs544719440   A   G    2
        576511      17  19159733  rs540831825   A   G    2
        717227      17  27196477  rs202111951   T   C   10
        919995      17  34642425  rs568794696   C   A    3
        2105598     17  77663493  rs148485780   C   T    5
                 CHROM       POS           ID REF ALT  VAR
        56637       17   1649616  rs544719440   A   C    1
        576510      17  19159733  rs540831825   A   C    1
        717226      17  27196477  rs202111951   T   A    9
        919587      17  34540858  rs568794696   C   A    3
        2105592     17  77663435  rs148485780   C   T    5        

        print(dataframe[dataframe.duplicated('ID', keep='first')])
        print(dataframe[dataframe.duplicated('ID', keep='last')])
 
        assert(len(dataframe['ID'].unique().tolist()) == len(dataframe['ID'].tolist()))
       
        """

        return dataframe

    def encode_ref_alt(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Returns: Dataframe which contains CHROM, POS, ID and VAR
        where VAR interprets the SNP variant given REF and ALT
        """

        def map(refalt) -> int:
            ref = refalt[0]
            alt = refalt[1]

            return DataHandler.mapping[ref][alt]

        dataframe['VAR'] = dataframe[['REF', 'ALT']].apply(map, axis=1)

        return dataframe

    def get_polyphen2_score(self, id: str) -> typing.Tuple[int, float]:
        """
        Returns: status of polyphen 2 score query and the corresponding polyphen 2 score for a given rsID
        """
        # TODO: get polyphen2 score
        return 0, 0.5

    def get_sift_score(self, id: str) -> typing.Tuple[int, float]:
        """
        Returns: status of sift score query and the corresponding sift score for a given rsID
        """
        # TODO: get sift score
        return 0, 0.5

    def get_cadd_score(self, id: str) -> typing.Tuple[int, float]:
        """
        Returns: status of CADD score query and the corresponding CADD score for a given rsID
        """
        # TODO: get cadd score
        return 0, 0.5

    def translate_chrom(self, id) -> int:
        """
        translate non integer chromosomes (X,Y & MT) to integers (23, 24 & 25)
        """
        if id['CHROM'] == 'X':
            return 23
        elif id['CHROM'] == 'Y':
            return 24
        elif id['CHROM'] == 'MT':
            return 25
        else:
            return int(id["CHROM"])


if __name__ == '__main__':
    URLS = {
        'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control.vcf",
        'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed.vcf",
    }

    from pathlib import Path

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()

    dh = DataHandler()

    trans_vcf = dh.translate_vcf("control_v2.vcf")

    trans_cvar = dh.translate_clinvar("clinvar.vcf")

    print(trans_vcf)
    print(trans_cvar)
