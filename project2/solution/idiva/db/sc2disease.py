# RA, 2020-11-23
# RA, 2020-12-05

from idiva import log

import pandas

import contextlib
import typing
import re
import pandas
import numpy
import idiva.io
import idiva.utils

from idiva import log
from tqdm import tqdm

URL = "http://easybioai.com/sc2disease/static/allgwas.txt"


def allgwas_reference() -> pandas.DataFrame:
    from idiva.download import download
    log.info(F"Downloading {URL}")
    with download(URL).now.open() as fd:
        log.info(F"Processing.")
        s: pandas.Series
        df = pandas.read_csv(fd, sep='\t', names=["gene", "rs", "p", "chrom", "pos", "disease"])
        df = pandas.DataFrame(
            columns=["ID", "SC2D"],
            data=[
                (rs, disease)
                for (rs, disease) in zip(df.rs, df.disease)
                for rs in re.findall(r"rs[0-9]+", rs)
            ]
        )
        df = df.groupby('ID', as_index=False).agg({'SC2D': lambda s: '"' + ", ".join(sorted(set(s))) + '"'})
        return df


def allgwas(*, case: idiva.io.ReadVCF, ctrl: idiva.io.ReadVCF):
    del ctrl

    ref = allgwas_reference()

    from idiva.clf.df import c5_df, join
    case = c5_df(case)

    class response:
        id_cols = list(case.columns)

        # One entry per columns, keyed by column name
        # https://samtools.github.io/hts-specs/VCFv4.1.pdf
        info = {
            'SC2D': {
                'Number': "1",
                'Type': "String",
                'Description': '"SC2Disease phenotype (http://easybioai.com/sc2disease/download)"'
            },
        }

        df = join(case=case, ctrl=ref, how='left', on='ID')

    catch: pandas.DataFrame = response.df[~response.df.SC2D.isna()]
    catch = catch.ID.groupby(catch.SC2D).agg(lambda s: ",".join(s))
    log.info(F"SC2Disease found: {dict(catch)}")

    return response


if __name__ == '__main__':
    from collections import Counter

    allgwas = allgwas_reference()

    # https://www.snpedia.com/index.php/Rs7903146
    print(allgwas[allgwas.ID == "rs7903146"])
