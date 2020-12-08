# LB 23-11-2020
import os
import typing
from pathlib import Path
import shlex
from subprocess import Popen

import numpy as np
import pandas as pd
from tqdm import tqdm

from idiva.fextr import FeatureExtractor, align
from idiva.io import ReadVCF
from idiva.utils import seek_then_rewind
from idiva import log
from idiva.clf.utils import TrainPhenomenetArgs
from idiva.db import db
import xlrd

MAPPING = {'A': {'C': 1, 'G': 2, 'T': 2},
           'C': {'A': 3, 'G': 4, 'T': 5},
           'G': {'A': 6, 'C': 7, 'T': 8},
           'T': {'A': 9, 'C': 10, 'G': 11}}


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

    mapping = MAPPING

    INIT_COLS = ["CHROM", "POS", "ID", "REF", "ALT"]
    CLINVAR_COLS = INIT_COLS
    CLINVAR_COLS_IDX = [0, 1, 2, 3, 4]

    def get_clf_datalines(self, df_clinvar: pd.DataFrame):
        """
        HK, 2020-11-21
        """
        for idx, row in tqdm(df_clinvar.iterrows(), total=len(df_clinvar), postfix='iterating df_clinvar'):
            if (str(row.ref) in ['A', 'C', 'G', 'T']) and (str(row.alt) in ['A', 'C', 'G', 'T']):
                line = {
                    'ID': row.id,
                    'CHROM': self.translate_chrom(row.chrom),
                    'POS': row['pos'],
                    'REF': row.ref,
                    'ALT': row.alt,
                    'VAR': self.mapping[row.ref][row.alt],
                    'label': 1 if row['CLNSIG'] == 'Pathogenic' else 0,

                }

                yield line

    def df_clinvar_to_clf_data(self, df_clinvar: pd.DataFrame) -> pd.DataFrame:
        """
        HK, 2020-11-21
        """
        dataframe = pd.DataFrame(data=self.get_clf_datalines(df_clinvar))
        dataframe = dataframe.drop_duplicates()
        dataframe = dataframe.sort_values(by=['CHROM', 'POS'])
        dataframe['CHROMPOSALTID'] = dataframe[['CHROM', 'POS', 'ALT']].apply(self.index_map, axis=1)
        dataframe = dataframe.set_index('CHROMPOSALTID')
        dataframe = dataframe.drop(columns=['REF', 'ALT'])

        return dataframe

    def get_phenomenet_training_data(self, args: TrainPhenomenetArgs) -> pd.DataFrame:
        if args.database == 'clinvar_dbSNP':
            clf_data = db.get_db_label_df(which_dbSNP=17, with_chrom_col=True).rename(columns={'class': 'label'})
            clf_data = clf_data.dropna(subset=['ref', 'alt'])
            clf_data = clf_data[clf_data.ref != 'N']
            clf_data = clf_data[clf_data.alt != 'N']
            clf_data['var'] = clf_data[['ref', 'alt']].apply(lambda x: self.mapping[x[0]][x[1]], axis=1)
        elif args.database == 'clinvar_processed':
            clf_data = self.get_clinvar_clf_processed_data()
        else:
            raise NotImplementedError(f"database {args.database} is not implemented.")

        return clf_data

    def get_clinvar_clf_data(self, clinvar_file: str = 'vcf_37') -> pd.DataFrame:
        """
        Loads clinvar_clf_data suitable for a classifier.
        Looks for in "_cache" or creates if not found two files:
            - the "exploded" clinvar file as a dataframe compatible csv (exploded meaning that all information from
            the INFO column is extracted to its own column)

            - the dataframe compatible csv file containing all extracted and encoded features
            to train a classifier from the clinvar dataframe

        HK, 2020-11-22
        RA, 2020-11-22
        """

        from idiva.io import cache_df

        def maker_clinvar() -> pd.DataFrame:
            from idiva.db import clinvar_open
            from idiva.io import ReadVCF
            from idiva.db.clinvar import clinvar_to_df

            with clinvar_open(which=clinvar_file) as fd:
                return clinvar_to_df(ReadVCF(fd))

        df_clinvar = cache_df(name=("clinvar_" + clinvar_file), key=[clinvar_file], df_maker=maker_clinvar)
        df_clinvar = df_clinvar.sort_values(by=['chrom', 'pos'])

        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]

        return cache_df(name="clinvar_clf_data", key=[clinvar_file, "v01"],
                        df_maker=lambda: self.df_clinvar_to_clf_data(df_clinvar_reduced))

    def get_clinvar_clf_processed_data(self) -> pd.DataFrame:

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        file_path = os.path.join(cache, 'training.csv')

        file_name = str(cache) + '/training.csv'

        print(os.path.isfile(file_path))

        if not os.path.isfile(file_path):
            self.preprocess_clinvar()

        dataframe = pd.read_csv(file_name, sep='\t', comment='#')

        return dataframe

    def get_clinvar_clf_extended(self, which: str = 'vcf_37'):
        from idiva.io import cache_df

        def maker_clinvar() -> pd.DataFrame:
            from idiva.db import clinvar_open
            from idiva.io import ReadVCF
            from idiva.db.clinvar import clinvar_to_df

            with clinvar_open(which=which) as fd:
                return clinvar_to_df(ReadVCF(fd))

        df_clinvar = cache_df(name=("clinvar_" + which), key=[which], df_maker=maker_clinvar)
        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]

        df_clinvar_reduced = df_clinvar_reduced.rename(columns={'chrom': 'CHROM', 'pos': 'POS',
                                                                'id': 'ID', 'ref': 'REF', 'alt': 'ALT', 'qual': 'QUAL',
                                                                'filter': 'FILTER'})

        # remove indels
        df_clinvar_reduced = df_clinvar_reduced[
            df_clinvar_reduced['REF'].apply(lambda x: str(x) in ['A', 'C', 'G', 'T'])]
        df_clinvar_reduced = df_clinvar_reduced[
            df_clinvar_reduced['ALT'].apply(lambda x: str(x) in ['A', 'C', 'G', 'T'])]

        df_clinvar_reduced = df_clinvar_reduced.drop_duplicates()

        df_clinvar_reduced['INFO'] = "."

        df_clinvar_reduced['labels'] = df_clinvar_reduced['CLNSIG'].apply(lambda row: 1 if row == 'Pathogenic' else 0)

        df_clinvar_reduced = df_clinvar_reduced[
            ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'labels']]

        return df_clinvar_reduced

    def create_training_set(self) -> typing.Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Returns training features and corresponding labels given a clinvar vcf file
        """

        clinvar_clf_data = self.get_clinvar_clf_processed_data()

        x_train = clinvar_clf_data.loc[:, clinvar_clf_data.columns != 'label']
        y_train = clinvar_clf_data.loc[:, clinvar_clf_data.columns == 'label']
        y_train.set_index(x_train.index)

        return x_train, y_train

    def create_test_set(self, vcf1, vcf2=None) -> pd.DataFrame:

        frame = self.translate_vcf(vcf1)

        if vcf2 is not None:
            baseframe2 = self.translate_vcf(vcf2)

            # merge both frames into one
            frame = pd.concat([frame, baseframe2])
            frame.drop_duplicates()

        return frame

    def create_test_set_v2(self, vcf_ctrl, vcf_case) -> pd.DataFrame:
        """
        creates test set by first reducing the number of samples and then adding sift and cadd scores
        """
        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        file_name = 'test.csv'

        file_path = os.path.join(cache, file_name)

        file_name = str(cache) + '/' + file_name

        dataframe_base = None

        # if file does not exists
        if not os.path.isfile(file_path):
            log.info("Annotate test set ")

            fextr = FeatureExtractor(ctrl_vcf=vcf_ctrl, case_vcf=vcf_case)
            dataframe_base = fextr.get_reduced_dataframe()

            dataframe_sift = self.add_sift_score(dataframe_base, 'our')

            dataframe_sift['CHROM'] = pd.to_numeric(dataframe_sift[['CHROM']].apply(self.translate_chrom, axis=1))

            dataframe_cadd = dataframe_base

            dataframe_cadd = dataframe_cadd[['CHROM', 'POS', 'ID', 'REF', 'ALT']]

            dataframe_cadd = self.add_cadd_score(dataframe_cadd)

            dataframe = dataframe_cadd
            dataframe[['SIFT_SCORE', 'SIFT_SUCC']] = dataframe_sift[['SIFT_SCORE', 'SIFT_SUCC']]

            dataframe['SIFT_SUCC'] = dataframe['SIFT_SUCC'].fillna(value=0)
            dataframe['SIFT_SCORE'] = dataframe['SIFT_SCORE'].fillna(value=0.05)
            dataframe['CADD_SUCC'] = dataframe['CADD_SUCC'].fillna(value=0)
            dataframe['CADD_PHRED'] = dataframe['CADD_PHRED'].fillna(value=30)

            dataframe = self.encode_ref_alt(dataframe)

            dataframe = dataframe[['CHROM', 'POS', 'VAR', 'CADD_PHRED', 'CADD_SUCC', 'SIFT_SCORE', 'SIFT_SUCC']]

            cols = ['CHROM', 'POS', 'VAR', 'CADD_PHRED', 'CADD_SUCC', 'SIFT_SCORE', 'SIFT_SUCC']

            dataframe[cols] = dataframe[cols].apply(pd.to_numeric, errors='coerce', axis=1)

            dataframe.to_csv(file_name, sep='\t')

        # load stored test set
        else:
            log.info("load stored test set")
            dataframe = pd.read_csv(file_name, sep='\t')

        return dataframe

    def translate_vcf(self, vcf) -> pd.DataFrame:
        """
        Returns a dataframe that contains the following features from a vcf file
        CHROM, POS, ID, VAR
        """

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        with ReadVCF.open(vcf) as reader:

            with seek_then_rewind(reader.fd, seek=reader.dataline_start_pos) as fd:

                dataframe = pd.read_csv(fd, sep='\t', usecols=range(len(DataHandler.INIT_COLS)), header=None,
                                        names=DataHandler.INIT_COLS,
                                        dtype={'CHROM': np.int, 'POS': np.int, 'ID': np.str, 'REF': np.str,
                                               'ALT': np.str})

                # Check if ALT contains only one value or several values seperated by ','
                assert (len([uni for uni in dataframe['ALT'].unique().tolist() if ',' in uni]) == 0)

                # store only SNP variants
                dataframe = dataframe[dataframe['REF'].apply(lambda x: {x}.issubset({'A', 'C', 'G', 'T'}))]
                dataframe = dataframe[dataframe['ALT'].apply(lambda x: {x}.issubset({'A', 'C', 'G', 'T'}))]

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
        #              consequence of real world data (Kjong Nov 30)
        #              => identify samples by CHROM, POS and VAR
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

    def add_sift_score(self, dataframe: pd.DataFrame, type: str) -> pd.DataFrame:
        """
        Appends sift scores and success to dataframe
        The dataframe needs at least following columns: CHROM, POS, REF, ALT
        """
        log.info("creating sift scores for " + type)
        # https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar
        # make dataframe compatible for sift annotator
        if 'ID' not in dataframe:
            dataframe['ID'] = '.'
        if 'QUAL' not in dataframe:
            dataframe['QUAL'] = '.'
        if 'FILTER' not in dataframe:
            dataframe['FILTER'] = '.'
        if 'INFO' not in dataframe:
            dataframe['INFO'] = '.'

        dataframe = dataframe.fillna(value=".")

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        file_path = str(cache) + "/" + type + "_sift.vcf"

        # create vcf file
        dataframe.rename(columns={'CHROM': '#CHROM'}).to_csv(file_path, sep='\t', index=False)

        file_name = type

        sift_folder = str(cache) + "/" + file_name + "_sift"

        cmd = 'java -jar ' + str(cache) + '/SIFT4G_Annotator.jar -c -i ' + file_path + ' -d ' + str(
            cache) + '/GRCh37.74 -r ' + sift_folder

        args = shlex.split(cmd)

        process = Popen(args)
        process.wait()

        sift_file = sift_folder + "/" + file_name + "_sift_SIFTannotations.xls"

        sift_dataframe = pd.read_table(sift_file)

        dataframe['SIFT_SCORE'] = np.nan
        dataframe['SIFT_SUCC'] = np.nan

        sift_iter = sift_dataframe.iterrows()
        next_sift = next(sift_iter, None)

        for idx, row in tqdm(dataframe.iterrows(), total=len(dataframe), postfix='inserting sift scores'):

            if next_sift is not None:
                if row['POS'] == next_sift[1]['POS']:
                    if not pd.isna(next_sift[1]['SIFT_SCORE']):
                        sift_score = next_sift[1]['SIFT_SCORE']
                        sift_succ = 1
                    else:
                        sift_score = np.nan
                        sift_succ = 0

                    next_sift = next(sift_iter, None)
                else:
                    sift_score = np.nan
                    sift_succ = 0
            else:
                sift_score = np.nan
                sift_succ = 0

            dataframe.loc[idx, 'SIFT_SCORE'] = sift_score
            dataframe.loc[idx, 'SIFT_SUCC'] = sift_succ

        return dataframe

    def add_cadd_score(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Returns: phred score and query status for a given dataframe containing CHROM, POS, REF, ALT
        """

        def fetch(chrom: int, poss: list, refs: list, alts: list, pbar: tqdm):

            from subprocess import Popen, PIPE

            next = {'A': {'start': 'C', 'C': 'G', 'G': 'T', 'T': 'start'},
                    'C': {'start': 'A', 'A': 'G', 'G': 'T', 'T': 'start'},
                    'G': {'start': 'A', 'A': 'C', 'C': 'T', 'T': 'start'},
                    'T': {'start': 'A', 'A': 'C', 'C': 'G', 'G': 'start'}}

            scores = np.empty(shape=(len(poss), 1))
            phreds = np.empty(shape=(len(poss), 1))
            succ = np.zeros(shape=(len(poss), 1))

            current_pos = 0

            old_ref = ''
            new_ref = ''

            current_alt = ''

            status = 0
            count = 0

            # in case of an "error" ([E::hts_open_format] Failed to open ...) retry (max 100 times)
            while status == 0 and count < 100:

                process = Popen(
                    ['tabix',
                     'https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs_inclAnno.tsv.gz',
                     'IndexFile', str(chrom) + ':' + str(poss[0]) + '-' + str(poss[-1])], stdout=PIPE)

                for idx, line in enumerate(process.stdout):

                    # Success, we got an output
                    status = 1

                    # interpret line as list
                    string_list = line.decode("utf-8").strip().split("\t")
                    # set new reference
                    new_ref = string_list[2]

                    # if we get the same information for the same alternative variant then skip it
                    if string_list[3] == current_alt:
                        continue
                    # set alternative to next nucleotide
                    else:
                        if old_ref != new_ref or next[new_ref][current_alt] == 'start':
                            current_alt = 'start'

                        current_alt = next[string_list[2]][current_alt]

                    # if the current line contains information about a asked position then store it
                    if int(string_list[1]) == poss[current_pos] and current_alt == alts[current_pos]:
                        """
                        # for debugging:
                        print("score", string_list[-2], "phred", string_list[-1], "pos", string_list[1], "ref", new_ref,
                              "alt", current_alt)
                        """

                        scores[current_pos] = string_list[-2]
                        phreds[current_pos] = string_list[-1]
                        succ[current_pos] = 1

                        current_pos += 1

                        # break if we iterated over all positions
                        if current_pos == len(poss):
                            break

                    old_ref = new_ref

                process.terminate()
                count += 1

            pbar.update(1)

            return [scores, phreds, succ]

        import concurrent.futures
        import multiprocessing

        chroms = dataframe['CHROM'].tolist()
        chroms = [self.translate_chrom_back(chrom) for chrom in chroms]

        poss = dataframe['POS'].tolist()
        refs = dataframe['REF'].tolist()
        alts = dataframe['ALT'].tolist()

        # start index of first bucket
        start_pos = 0
        buckets = [start_pos]

        current_chrom = chroms[0]

        # bucket spans a range of max cut_off nucleotide bases
        cut_off = 50000

        for idx, pos in enumerate(poss[1:], 1):
            if pos - poss[start_pos] > cut_off or current_chrom != chroms[idx]:
                start_pos = idx
                current_chrom = chroms[idx]
                buckets.append(idx)

        buckets.append(len(poss))

        futures = []

        with tqdm(total=len(buckets) - 1, postfix='creating CADD scores') as pbar:
            with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count() - 1) as executor:
                for idx in range(len(buckets) - 1):
                    args = [chroms[idx],
                            poss[buckets[idx]:buckets[idx + 1]],
                            refs[buckets[idx]:buckets[idx + 1]],
                            alts[buckets[idx]:buckets[idx + 1]],
                            pbar]

                    futures.append(executor.submit(lambda p: fetch(*p), args))

        scores = futures[0].result()[0].ravel()
        phreds = futures[0].result()[1].ravel()
        succ = futures[0].result()[2].ravel()

        for idx, future in enumerate(futures[1:], 1):
            scores = np.concatenate([scores, future.result()[0].ravel()])
            phreds = np.concatenate([phreds, future.result()[1].ravel()])
            succ = np.concatenate([succ, future.result()[2].ravel()])

        dataframe['CADD_PHRED'] = phreds
        dataframe['CADD_SUCC'] = succ

        return dataframe

    def translate_chrom(self, chrom: typing.Union[str, int]) -> int:
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

    def translate_chrom_back(self, chrom: int) -> typing.Union[str, int]:
        if chrom == 23:
            return 'X'
        elif chrom == 24:
            return 'Y'
        elif chrom == 25:
            return 'MT'
        else:
            return chrom

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

    def preprocess_clinvar(self, which='vcf_37'):

        log.info("preprocessing clinvar file")

        # TODO:
        #  How to download clinvar file, unzip it and store it as clinvar.vcf in download_cache???

        # TODO:
        #  How to download GRCH37.74 unzip it and store it in download_cache???

        dataframe_base = self.get_clinvar_clf_extended()
        dataframe_base = dataframe_base.reset_index(drop=True)

        labels = dataframe_base['labels']
        dataframe_base = dataframe_base.drop('labels', axis=1)

        dataframe_sift = self.add_sift_score(dataframe_base, 'clinvar')

        dataframe_sift['CHROM'] = pd.to_numeric(dataframe_sift[['CHROM']].apply(self.translate_chrom, axis=1))

        dataframe_cadd = dataframe_base

        dataframe_cadd['CHROM'] = pd.to_numeric(dataframe_cadd[['CHROM']].apply(self.translate_chrom, axis=1))
        dataframe_cadd = dataframe_cadd[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
        dataframe_cadd = dataframe_cadd.reset_index(drop=True)

        cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
        assert cache.is_dir()

        # get cadd annotations
        from idiva.download import download

        data = download('https://polybox.ethz.ch/index.php/s/GRgDYHOAaw75D60/download').now
        with data.open() as fd:
            cadd_scores = pd.read_csv(fd, sep='\t', usecols=range(1, 6), comment='#',
                                      names=['CHROM', 'POS', 'REF', 'ALT', 'CADD_SCORE', 'CADD_PHRED'])

        cadd_scores['CADD_SUCC'] = 1

        dataframe_cadd[['CADD_PHRED', 'CADD_SUCC']] = cadd_scores[['CADD_PHRED', 'CADD_SUCC']]

        dataframe = dataframe_cadd
        dataframe[['SIFT_SCORE', 'SIFT_SUCC']] = dataframe_sift[['SIFT_SCORE', 'SIFT_SUCC']]

        dataframe['SIFT_SCORE'] = dataframe['SIFT_SCORE'].fillna(value=0.05)
        dataframe['CADD_SUCC'] = dataframe['CADD_SUCC'].fillna(value=0)
        dataframe['CADD_PHRED'] = dataframe['CADD_PHRED'].fillna(value=30)

        dataframe = self.encode_ref_alt(dataframe)

        dataframe = dataframe[['CHROM', 'POS', 'VAR', 'CADD_PHRED', 'CADD_SUCC', 'SIFT_SCORE', 'SIFT_SUCC']]

        dataframe.insert(7, 'label', labels)

        cols = ['CHROM', 'POS', 'VAR', 'CADD_PHRED', 'CADD_SUCC', 'SIFT_SCORE', 'SIFT_SUCC', 'label']

        dataframe[cols] = dataframe[cols].apply(pd.to_numeric, errors='coerce', axis=1)

        file_path = str(cache) + "/training.csv"

        dataframe.to_csv(file_path, sep='\t', index=False)

        return dataframe


if __name__ == '__main__':
    dh = DataHandler()

    # print(dh.preprocess_clinvar())

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()

    with ReadVCF.open(str(cache) + '/control_v2.vcf') as ctrl_vcf:
        with ReadVCF.open(str(cache) + '/case_processed_v2.vcf') as case_vcf:
            test_set = dh.create_test_set_v2(ctrl_vcf, case_vcf)


    """
    print(dataframe)

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()

    file_path = str(cache) + "/cadd_full.vcf"

    dataframe = dataframe.fillna(value=".")

    dataframe.rename(columns={'CHROM': '#CHROM'}).to_csv(file_path, sep='\t', index=False)
    """
