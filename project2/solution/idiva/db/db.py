# HK, 2020-12-05

import pandas as pd
from idiva import log


def join_clinvar_dbSNP(df_clinvar: pd.DataFrame, df_dbSNP: pd.DataFrame) -> pd.DataFrame:
    merge_on_pos_ref_alt = (
        df_clinvar[['pos', 'ref', 'alt', 'CLNSIG']]
            .replace({'Pathogenic': 1, 'Benign': 0})
            .merge(df_dbSNP[['pos', 'ref', 'alt', 'CLNSIG']]
                   .replace({'CLNSIG': {2: 0, 5: 1}}),
                   left_on=['pos', 'ref', 'alt'],
                   right_on=['pos', 'ref', 'alt'],
                   how='outer',
                   suffixes=('_clinvar', '_dbSNP'),
                   )
    )

    merge_on_pos_ref_alt['class'] = merge_on_pos_ref_alt['CLNSIG_clinvar'].fillna(0).astype(int) & merge_on_pos_ref_alt[
        'CLNSIG_dbSNP'].fillna(0).astype(int)
    merge_on_pos_ref_alt = merge_on_pos_ref_alt.drop(['CLNSIG_clinvar', 'CLNSIG_dbSNP'], axis=1)

    return merge_on_pos_ref_alt


def get_db_label_df(clinvar_file: str = 'vcf_37', which_dbSNP: int = 17) -> pd.DataFrame:
    """
    which_dbSNP: integer indicating which chromosome is looked up in the dbSNP
    Returns a dataframe like the following:

                       pos ref alt  class

    0           949422   G   A    0.0
    1           949422   G   A    0.0
    2           949438   A   G    0.0
    3           949438   A   G    0.0
    4           949523   C   T    1.0

    The class label is obtained by joining the clinvar dataframe with the dbSNP dataframe on pos, ref, alt.
    the class label from left and right are joined by taking a logical AND of both.
    """
    from idiva.db.dbSNP import get_dbSNP_df
    from idiva.io import cache_df
    log.info("Getting database labels")

    def maker_merge_on_pos_ref_alt() -> pd.DataFrame:
        log.info("Creating database labels dataframe")
        dbSNP_df = get_dbSNP_df(which_dbSNP)
        reduced_dbSNP = dbSNP_df.loc[(dbSNP_df.CLNSIG == 2) | (dbSNP_df.CLNSIG == 5)]

        def maker_clinvar() -> pd.DataFrame:
            """
            creates the clinvar dataframe
            """
            from idiva.db import clinvar_open
            from idiva.io import ReadVCF
            from idiva.db.clinvar import clinvar_to_df

            with clinvar_open(which=clinvar_file) as fd:
                return clinvar_to_df(ReadVCF(fd))

        df_clinvar = cache_df(name=("clinvar_" + clinvar_file), key=[clinvar_file], df_maker=maker_clinvar)
        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]

        merge_on_pos_ref_alt = join_clinvar_dbSNP(df_clinvar=df_clinvar_reduced, df_dbSNP=reduced_dbSNP)

        return merge_on_pos_ref_alt.astype({'pos': int, 'ref': 'category', 'alt': 'category', 'class': int})

    merge_on_pos_ref_alt = cache_df(name=("db_PosRefAlt" + clinvar_file), key=[clinvar_file],
                                    df_maker=maker_merge_on_pos_ref_alt)

    return merge_on_pos_ref_alt
