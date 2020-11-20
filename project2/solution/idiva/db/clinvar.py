# RA, 2020-11-11

import contextlib
import gzip
import io
import typing
from itertools import product

import pandas as pd
from tqdm import tqdm

import idiva.utils
from idiva.clf.utils import NucEncoder
from idiva.io.vcf import ReadVCF
from idiva.utils import at_most_n

URL = {
    'vcf_37': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    'vcf_38': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}

# pandas cannot represent nan as integers.
dtype_clinvar_df = {'chrom': str, 'pos': pd.Int64Dtype(), 'id': pd.Int64Dtype(), 'ref': str, 'alt': str, 'qual': str,
                    'filter': str,
                    'format': str, 'samples': str, 'ALLELEID': pd.Int64Dtype(), 'CLNDISDB': str, 'CLNDN': str,
                    'CLNHGVS': str,
                    'CLNREVSTAT': str, 'CLNSIG': str, 'CLNVC': str, 'CLNVCSO': str, 'GENEINFO': str,
                    'MC': str, 'ORIGIN': pd.Int64Dtype(), 'AF_ESP': float, 'AF_EXAC': float, 'AF_TGP': float,
                    'RS': pd.Int64Dtype()}


@contextlib.contextmanager
def clinvar_open(which='vcf_37') -> typing.Iterable[typing.TextIO]:
    from idiva.download import download
    data = download(URL[which]).now
    with data.open(mode='rb') as gz:
        with gzip.open(gz) as fd:
            yield io.TextIOWrapper(fd)


def clinvar_meta(which='vcf_37') -> idiva.utils.minidict:
    from idiva.download import download
    data = download(URL[which]).now
    return idiva.utils.minidict(data.meta)


def get_info_dict(info: str) -> dict:
    """
    Yields info dict for every RS id.
    The info field can contain several RS ids and several OMIM ids. This function yields a dict for any
    combination of them.
    """
    import re
    OMIM_ids = [None]
    RS_ids = [None]
    info_dict = {}
    # go through info which is semicolon separated
    for elem in info.split(';'):
        # spit into key, value
        k, v = elem.split('=')
        if k == 'CLNDISDB':
            OMIM_ids = re.findall('OMIM:\d+', v)
        info_dict[k] = v
        if k == 'RS':
            for rs_id in v.split('|'):
                RS_ids.append('rs' + str(int(rs_id)))

    for OMIM_id, RS_id in product(OMIM_ids, RS_ids):
        info_dict['OMIM_id'] = OMIM_id
        info_dict['RS'] = RS_id
        yield info_dict


def clf_datalines(df_clinvar: pd.DataFrame, nuc_encoder: NucEncoder):
    for idx, row in tqdm(df_clinvar.iterrows(), total=len(df_clinvar), postfix='iterating df_clinvar'):
        line = {
            'pos': row['pos'],
            'ref': nuc_encoder.encode(None) if str(row['alt']) == 'nan' else nuc_encoder.encode(row['ref']),
            'alt': nuc_encoder.encode(None) if str(row['alt']) == 'nan' else nuc_encoder.encode(row['alt']),
            'label': 1 if row['CLNSIG'] == 'Pathogenic' else 0
        }
        yield line


def df_clinvar_to_clf_data(df_clinvar: pd.DataFrame):
    nuc_encoder = NucEncoder()
    return pd.DataFrame(data=clf_datalines(df_clinvar, nuc_encoder))


def clinvar_datalines(vcf: idiva.io.ReadVCF):
    for idx, line in tqdm(enumerate(vcf.datalines), postfix='reading clinvar file'):
        for info_dict in get_info_dict(line.info):
            line_dict = {k: line.__dict__[k] for k in line.__dict__.keys() if not k == 'info'}
            line_dict = dict(line_dict, **info_dict)

            yield line_dict


def clinvar_to_df(vcf: idiva.io.ReadVCF) -> pd.DataFrame:
    """
    Creates a dataframe from the clinvar file. Adds all the INFO fields as additional columns.
    """

    return pd.DataFrame(data=clinvar_datalines(vcf)).astype({'ref': str, 'alt': str})


def clinvar_rs_ids(which='vcf_37'):
    from idiva.io.vcf import ReadVCF
    with clinvar_open(which) as fd:
        for dataline in ReadVCF(fd):
            pass


if __name__ == '__main__':
    with clinvar_open('vcf_37') as fd:
        reader = ReadVCF(fd)
        print(*at_most_n(reader.datalines, n=10), sep='\n')
