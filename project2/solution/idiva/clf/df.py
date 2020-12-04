from idiva import log

import idiva.io

import typing
import pandas
import numpy

# Do not remove:
import joblib
import tqdm


def apply_dtype(df: pandas.DataFrame) -> pandas.DataFrame:
    dtype_v0 = {'CHROM': str, 'POS': int, 'ID': str, 'ALT0': float, 'ALT1': float, 'ALT2': float}
    return df.astype({c: dtype_v0[c] for c in df.columns if c in dtype_v0})


def join(*, case: pandas.DataFrame, ctrl: pandas.DataFrame, how="outer", on=['CHROM', 'POS', 'ID']) -> pandas.DataFrame:
    """
    Outer-join two dataframes on the columns CHROM, POS, ID.
    Use the suffixes _case and _ctrl for the other ambiguous columns.

    RA, 2020-11-14
    LB, 2020-12-04
    RA, 2020-12-04
    """

    df = pandas.merge_ordered(
        left=case, right=ctrl,
        suffixes=['_case', '_ctrl'],
        on=on,
        how=how,
    )

    return df


def v0_datalines(vcf: idiva.io.ReadVCF) -> typing.Iterable[dict]:
    """
    RA, 2020-11-14
    RA, 2020-12-01
    """
    from idiva.io.vcf import parse_gt

    # special = sorted({F"<{k}>" for k in vcf.meta['ALT'].keys()})

    def parse(dataline):
        samples = tuple(parse_gt(gt) for gt in dataline.samples)

        # SIMPLIFICATION

        samples = tuple(((a != 0) + (b != 0)) for (a, b) in samples)
        # samples = (numpy.array(samples) != 0).sum(axis=1)

        line = {
            'CHROM': dataline.chrom,
            'POS': dataline.pos,
            'ID': dataline.id,
            'ALT0': sum((x == 0) for x in samples),
            'ALT1': sum((x == 1) for x in samples),
            'ALT2': sum((x == 2) for x in samples),
        }

        return line

    log.debug("Parsing datalines.")

    from tqdm import tqdm
    from joblib import Parallel, delayed

    lines = Parallel(n_jobs=8, prefer="processes")(
        delayed(parse)(dataline)
        for dataline in tqdm(vcf)
    )

    yield from lines


def c3_df(vcf: idiva.io.ReadVCF) -> pandas.DataFrame:
    """
    Returns a dataframe like

           CHROM   POS           ID
        0     17    52  rs556541063
        1     17    56  rs145615430
        2     17    78  rs148170422
        ..   ...   ...          ...

    Rewinds the `vcf`.

    RA, 2020-12-03
    """

    with vcf.rewind_when_done:
        return apply_dtype(pandas.DataFrame(
            data=((dataline.chrom, dataline.pos, dataline.id) for dataline in vcf),
            columns=["CHROM", "POS", "ID"],
        ))


def v0_df(vcf: idiva.io.ReadVCF) -> pandas.DataFrame:
    """
    Make a dataframe from the VCF datalines.
    Only counts the presence of alternative variants.
    Returns a dataframe like

           CHROM   POS           ID   ALT0  ALT1  ALT2
        0     17    52  rs556541063  314.0   0.0   0.0
        1     17    56  rs145615430  314.0   0.0   0.0
        2     17    78  rs148170422  313.0   1.0   0.0
        ..   ...   ...          ...    ...   ...   ...

    ID is not guaranteed to be unique.

    `int` data type is used for the ALT columns (see `dtype_v0`).

    Rewinds the `vcf`.

    RA, 2020-11-14
    """

    with vcf.rewind_when_done:
        log.debug("Getting the datalines and creating the dataframe.")
        return apply_dtype(pandas.DataFrame(data=v0_datalines(vcf)))


def get_clinvar_clf_data(base_string_encoding: str = 'integer') -> pandas.DataFrame:
    """
    Loads clinvar_clf_data suitable for a classifier.

    HK, 2020-11-22
    RA, 2020-11-22
    """

    from idiva.io import cache_df

    which = 'vcf_37'

    def maker_clinvar() -> pandas.DataFrame:
        from idiva.db import clinvar_open
        from idiva.io import ReadVCF
        from idiva.db.clinvar import clinvar_to_df

        with clinvar_open(which=which) as fd:
            return clinvar_to_df(ReadVCF(fd))

    df_clinvar = cache_df(name=("clinvar_" + which), key=[], df_maker=maker_clinvar)

    def maker_clinvar_clf() -> pandas.DataFrame:
        from idiva.db.clinvar import df_clinvar_to_clf_data
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
    from pathlib import Path
    from contexttimer import Timer

    cache = (Path(__file__).parent.parent.parent.parent / "input/download_cache").resolve()
    assert cache.is_dir()
    download = download.to(abs_path=cache)

    with download(URLS['ctrl']).now.open() as fd:
        with Timer() as timer:
            df = v0_df(idiva.io.ReadVCF(fd))

        print(F"This took {timer.elapsed} seconds")
