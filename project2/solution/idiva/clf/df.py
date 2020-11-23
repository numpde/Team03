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

    RA, 2020-11-14
    """

    df = pandas.merge_ordered(
        left=case, right=ctrl,
        suffixes=['_case', '_ctrl'],
        on=['CHROM', 'POS', 'ID'],
        how="outer",
    )

    return df


def v0_datalines(vcf: idiva.io.ReadVCF) -> typing.Iterable[dict]:
    """
    RA, 2020-11-14
    """
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

    Pandas Int64 data type is used for the ALT columns.

    RA, 2020-11-14
    """

    return pandas.DataFrame(data=v0_datalines(vcf)).astype(dtype_v0)


def get_clinvar_clf_data(base_string_encoding: str = 'integer') -> pd.DataFrame:
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

    def maker_clinvar_clf() -> pd.DataFrame:
        # If you change this function, change the cache key also.
        # Preparing the clinvar dataframe for categorical classification:
        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]
        return df_clinvar_to_clf_data(df_clinvar_reduced, base_string_encoding=base_string_encoding)

    return cache_df(name="clinvar_clf_data", key=[base_string_encoding, "v01"], df_maker=maker_clinvar_clf)


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
