# RA, 2020-12-01

import idiva.io


def vcf_to_fisher(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):
    from idiva import log

    from idiva.io import cache_df
    from idiva.clf.df import v0_df, join
    from idiva.stat import v0_fisher

    log.info("Creating a slim case/ctrl dataframe.")

    def df_maker1():
        # Note: v0_df is parallelized
        return join(case=v0_df(case), ctrl=v0_df(ctrl))

    cached_file_prefix = "case_ctrl__v0df"

    df = cache_df(cached_file_prefix, key=[case.md5, ctrl.md5], df_maker=df_maker1)
    log.debug(df)

    log.info("Computing the p-values case vs control.")

    def df_maker2():
        return v0_fisher(df).assign(CHROM=df.CHROM, POS=df.POS, ID=df.ID)

    pvalues = cache_df(cached_file_prefix + "__pvalues", key=[case.md5, ctrl.md5], df_maker=df_maker2)
    log.debug(pvalues)

    from idiva.clf.df import apply_dtype
    return apply_dtype(pvalues)
