# RA, 2020-12-01

import idiva.io
from idiva import log


def vcf_to_fisher(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):

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

    from idiva.clf.df import apply_dtype
    pvalues = apply_dtype(pvalues)

    pvalues.rename(inplace=True, columns={
        'ALT0_vs_Other': "F0",
        'ALT1_vs_Other': "F1",
        'ALT2_vs_Other': "F2",
    })

    log.debug(pvalues)

    class response:
        # One entry per columns, keyed by column name
        # https://samtools.github.io/hts-specs/VCFv4.1.pdf
        info = {
            'F0': {'Number': "1", 'Type': "Float", 'Description': '"Fisher test case vs ctrl on ALT0 vs rest"'},
            'F1': {'Number': "1", 'Type': "Float", 'Description': '"Fisher test case vs ctrl on ALT1 vs rest"'},
            'F2': {'Number': "1", 'Type': "Float", 'Description': '"Fisher test case vs ctrl on ALT2 vs rest"'},
        }

        df = pvalues

    return response



