# RA, 2020-12-01

import typing
import re
import pandas
import idiva.io
import idiva.utils

from idiva import log
from tqdm import tqdm


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
        id_cols = ["CHROM", "POS", "ID"]

        # One entry per columns, keyed by column name
        # https://samtools.github.io/hts-specs/VCFv4.1.pdf
        info = {
            'F0': {'Number': "1", 'Type': "Float", 'Description': '"Fisher test  case vs ctrl  x  0/0 vs rest"'},
            'F1': {'Number': "1", 'Type': "Float",
                   'Description': '"Fisher test  case vs ctrl  x  (0/1 and 1/0) vs rest"'},
            'F2': {'Number': "1", 'Type': "Float", 'Description': '"Fisher test  case vs ctrl  x  1/1 vs rest"'},
        }

        df = pvalues

    return response


def extract_pvalues(vcf: idiva.io.ReadVCF) -> pandas.DataFrame:
    from idiva.utils import unlist1

    with vcf.rewind_when_done:
        def data():
            log.info("Obtaining SC2Disease.")
            from idiva.db.sc2disease import allgwas_reference
            sc2d = dict(allgwas_reference()[["ID", "SC2D"]].set_index("ID").SC2D)

            log.info("Obtaining ClinVar.")
            from idiva.db.clinvar import clinvar_rsid_clnsig
            clnsig = dict(clinvar_rsid_clnsig()[["ID", "ClnSig"]].set_index("ID").ClnSig)

            log.info("Reading p-values from VCF.")
            for dataline in tqdm(vcf):
                try:
                    p = [
                        float(unlist1(re.findall(rF"{FX}=([^;]+);", dataline.info.strip() + ";")))
                        for FX in ["F0", "F1", "F2"]
                    ]
                    yield (
                        dataline.chrom, dataline.pos, dataline.id,
                        *p,
                        sc2d.get(dataline.id, "N/A"),
                        clnsig.get(dataline.id, "N/A"),
                    )
                except ValueError:
                    # Maybe encountered a '.'
                    pass

        return pandas.DataFrame(data=data(), columns=["CHROM", "POS", "ID", "F0", "F1", "F2", "SC2Disease", "ClnSig"])


def figure_pvalues(vcf: idiva.io.ReadVCF) -> typing.Iterable[idiva.utils.Plox]:
    import numpy as np
    import pandas as pd
    from idiva.utils import Plox

    pscores = extract_pvalues(vcf)

    for FX in ["F0", "F1", "F2"]:

        # with Plox() as px:
        #     px.a.hist(pscores[pscores != 1.0].values)
        #     yield px

        # pscores = pscores.sort_values(ascending=True).nsmallest(n=10000)

        # with Plox() as px:
        #     rank = (1 + np.arange(0, len(pscores)))
        #     px.a.plot(-np.log10(rank), -np.log10(pscores), '.-')
        #     # px.a.plot(-np.log10(expect), -np.log10(expect), '--', lw=0.5)
        #     px.a.grid()
        #     px.a.set_xlabel(r"$-\log_{10}$ rank")
        #     px.a.set_ylabel(r"$-\log_{10}$ p-value")
        #     yield px

        n = 2000
        pscores_fx = pscores.loc[pscores[FX].nsmallest(n=n).index]
        assert isinstance(pscores_fx, pd.DataFrame)

        for (chrom, df) in pscores_fx.groupby(pscores_fx.CHROM):
            with Plox() as px:
                spec = {'F0': "0/0 vs rest", 'F1': "(0/1 and 1/0) vs rest", 'F2': "1/1 vs rest"}[FX]
                px.info = {'name proposal': F"pvalues_chrom={chrom}_fx={FX}_n={n}", 'df': df, }
                px.a.plot(pscores_fx.POS, -np.log10(pscores_fx[FX]), '.')
                ylim = [0, max(max(px.a.get_ylim()), 1)]
                px.a.set_yticks(list(range(0, 2 + round(max(ylim)))))
                px.a.set_ylim(*ylim)
                px.a.set_xlabel(F"position on chr. {chrom}")
                px.a.set_ylabel(r"$-\log_{10}$ p-value")
                px.a.set_title(F"{n} smallest p-values: {spec}")
                px.a.grid()
                yield px
