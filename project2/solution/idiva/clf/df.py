# RA, 2020-11-14

import typing
from pathlib import Path

import pandas
import pandas as pd

import numpy
import idiva.io
from idiva.db.clinvar import df_clinvar_to_clf_data
from joblib import Parallel, delayed

dtype_v0 = {'CHROM': str, 'POS': int, 'ID': str, 'ALT0': float, 'ALT1': float, 'ALT2': float}


def join(*, case: pandas.DataFrame, ctrl: pandas.DataFrame) -> pandas.DataFrame:
    """
    Outer-join two dataframes on the columns CHROM, POS, ID.
    Use the suffixes _case and _ctrl for the other ambiguous columns.
    Uses sort=False in pandas.merge_ordered.
    """

    df = pandas.merge_ordered(
        left=case, right=ctrl,
        suffixes=['_case', '_ctrl'],
        on=['CHROM', 'POS', 'ID'],
        how="outer",
    )

    return df



def v0_datalines(vcf: idiva.io.ReadVCF) -> typing.Iterable[dict]:
    from idiva.io.vcf import parse_gt

    # special = sorted({F"<{k}>" for k in vcf.meta['ALT'].keys()})

    for dataline in vcf:
        samples = [parse_gt(gt) for gt in dataline.samples]

        # SIMPLIFICATION

        samples = (numpy.array(samples, dtype=int) != 0).sum(axis=1)

        line = {
            'CHROM': dataline.chrom,
            'POS': dataline.pos,
            'ID': dataline.id,
            'ALT0': (samples == 0).sum(),
            'ALT1': (samples == 1).sum(),
            'ALT2': (samples == 2).sum(),
        }

        yield line


def v0_df(vcf: idiva.io.ReadVCF) -> pandas.DataFrame:
    """
    Make a dataframe from the VCF datalines.
    Only counts the presence of alternative variants.
    Returns a dataframe like

        ID          CHROM   POS  ALT0  ALT1  ALT2
        rs556541063    17    52   500     0     0
        ...           ...   ...   ...   ...   ...
        rs191439637    17  2833   495     5     0
        rs149757548    17  2837   498     2     0

    ID is not guaranteed to be unique.

    Int64 data type is used for the ALT columns.
    """

    return pandas.DataFrame(data=v0_datalines(vcf)).astype(dtype_v0)


def get_clinvar_clf_data(data_dir: Path, save_df=False, base_string_encoding: str = 'integer') -> pd.DataFrame:
    """
    Loads clinvar_clf_data if it exists as csv under data_dir.
    Otherwise tries to load df_clinvar, and creates clinvar_clf_data from it.
    """
    # dict to store an index of the base_string_encoding for the file name endings
    base_string_encoding_idx = {'integer': 1, 'base_string_length': 2}

    clinvar_csv_path = data_dir / 'clinvar.csv.gzip'
    clinvar_clf_data_path = data_dir / f'clinvar_clf_data_{base_string_encoding_idx[base_string_encoding]}.csv.gzip'

    if not clinvar_clf_data_path.is_file():
        if not clinvar_csv_path.is_file():
            from idiva.db import clinvar_open
            from idiva.io import ReadVCF
            from idiva.db.clinvar import clinvar_to_df

            with clinvar_open(which='vcf_37') as fd:
                df_clinvar = clinvar_to_df(ReadVCF(fd))
            if save_df:
                df_clinvar.to_csv(clinvar_csv_path, index=False, compression='gzip')
        else:
            df_clinvar = pd.read_csv(clinvar_csv_path, compression='gzip')

        df_clinvar = df_clinvar.loc[(df_clinvar['CLNSIG'] == 'Pathogenic') | (df_clinvar['CLNSIG'] == 'Benign')]
        # transforming the clinvar dataframe to categorical data for the classifier
        clinvar_clf_data = df_clinvar_to_clf_data(df_clinvar, base_string_encoding=base_string_encoding)

        if save_df:
            clinvar_clf_data.to_csv(clinvar_clf_data_path, index=False, compression='gzip')
    else:
        clinvar_clf_data = pd.read_csv(clinvar_clf_data_path, compression='gzip')

    return clinvar_clf_data


if __name__ == '__main__':
    URLS = {
        'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
        'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
    }

    from tcga.utils import download
    from idiva.utils import at_most_n
    from pathlib import Path
    from contexttimer import Timer

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()
    download = download.to(abs_path=cache)

    with download(URLS['ctrl']).now.open() as fd:
        with Timer() as timer:
            df = v0_df(idiva.io.ReadVCF(fd))

        print(F"This took {timer.elapsed} seconds")
