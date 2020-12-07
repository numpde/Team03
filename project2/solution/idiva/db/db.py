# HK, 2020-12-05

import pandas as pd
from idiva import log

chromid_to_chromidx = {'NC_000001': 1, 'NC_000002': 2, 'NC_000003': 3, 'NC_000004': 4, 'NC_000005': 5, 'NC_000006': 6,
                       'NC_000007': 7, 'NC_000008': 8, 'NC_000009': 9, 'NC_000010': 10, 'NC_000011': 11,
                       'NC_000012': 12, 'NC_000013': 13, 'NC_000014': 14, 'NC_000015': 15, 'NC_000016': 16,
                       'NC_000017': 17, 'NC_000018': 18, 'NC_000019': 19, 'NC_000020': 20, 'NC_000021': 21,
                       'NC_000022': 22, 'NC_000023': 23, 'NC_000024': 24}


def join_clinvar_dbSNP(df_clinvar: pd.DataFrame, df_dbSNP: pd.DataFrame, with_chrom_col: bool = False) -> pd.DataFrame:
    cols = ['pos', 'ref', 'alt', 'CLNSIG']
    on = ['pos', 'ref', 'alt', 'chrom']
    if with_chrom_col:
        for k, v in chromid_to_chromidx.items():
            df_dbSNP['chrom'] = df_dbSNP['chrom'].replace(to_replace=f'^{k}', value=v, regex=True)
        cols.append('chrom')
        on.append('chrom')
        # todo are X,Y the 23rd and 24th chrom in dbSNP? Could be better to remove them if not sure.
        #  Also is MT chrom needed?
        df_clinvar['chrom'] = df_clinvar['chrom'].replace({'X': 23, 'Y': 24, 'MT': 25})

    merge_on_pos_ref_alt = (
        df_clinvar[cols].replace({'Pathogenic': 1, 'Benign': 0})
            .merge(df_dbSNP[cols]
                   .replace({'CLNSIG': {2: 0, 5: 1}}),
                   on=on,
                   how='outer',
                   suffixes=('_clinvar', '_dbSNP'),
                   )
    )

    merge_on_pos_ref_alt['class'] = merge_on_pos_ref_alt['CLNSIG_clinvar'].fillna(0).astype(int) & merge_on_pos_ref_alt[
        'CLNSIG_dbSNP'].fillna(0).astype(int)
    merge_on_pos_ref_alt = merge_on_pos_ref_alt.drop(['CLNSIG_clinvar', 'CLNSIG_dbSNP'], axis=1)

    return merge_on_pos_ref_alt


def get_db_label_df(clinvar_file: str = 'vcf_37', which_dbSNP: int = 17, with_chrom_col: bool = False) -> pd.DataFrame:
    """
    which_dbSNP: integer indicating which chromosome is looked up in the dbSNP
    with_chrom_col: boolean indicating if the returned dataframe should contain the chrom column
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
            log.info('Making clinvar df.')
            with clinvar_open(which=clinvar_file) as fd:
                return clinvar_to_df(ReadVCF(fd))

        df_clinvar = cache_df(name=("clinvar_" + clinvar_file), key=[clinvar_file], df_maker=maker_clinvar)
        df_clinvar_reduced = df_clinvar[df_clinvar['CLNSIG'].isin({'Pathogenic', 'Benign'})]

        merge_on_pos_ref_alt = join_clinvar_dbSNP(df_clinvar=df_clinvar_reduced, df_dbSNP=reduced_dbSNP,
                                                  with_chrom_col=with_chrom_col)

        types = {'pos': int, 'ref': 'category', 'alt': 'category', 'class': int}
        if with_chrom_col:
            types['chrom'] = int
        return merge_on_pos_ref_alt.astype(types)

    from idiva.clf.df import apply_dtype
    merge_on_pos_ref_alt = apply_dtype(
        cache_df(
            name=("db_PosRefAlt" + clinvar_file),
            key=[clinvar_file, str(which_dbSNP), str(with_chrom_col)],
            df_maker=maker_merge_on_pos_ref_alt
        )
    )

    return merge_on_pos_ref_alt
