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
    """
    RA, 2020-11-11
    """
    from idiva.download import download
    data = download(URL[which]).now
    with data.open(mode='rb') as gz:
        with gzip.open(gz) as fd:
            yield io.TextIOWrapper(fd)


def clinvar_meta(which='vcf_37') -> idiva.utils.minidict:
    """
    RA, 2020-11-11
    """
    from idiva.download import download
    data = download(URL[which]).now
    return idiva.utils.minidict(data.meta)


def get_info_dict(info: str) -> dict:
    """
    Yields info dict for every RS id.
    The info field can contain several RS ids and several OMIM ids. This function yields a dict for any
    combination of them.

    HK, 2020-11-12
    """
    import re
    OMIM_ids = []
    RS_ids = []
    info_dict = {}
    # go through info which is semicolon separated
    for elem in info.split(';'):
        # spit into key, value
        k, v = elem.split('=')
        if k == 'CLNDISDB':
            OMIM_ids = re.findall(r'OMIM:\d+', v)
        info_dict[k] = v
        if k == 'RS':
            for rs_id in v.split('|'):
                RS_ids.append('rs' + str(int(rs_id)))

    yield info_dict
    """
    # LB 6-12-2020 adapted
    if len(OMIM_ids) > 0 and len(RS_ids) > 0:
        for OMIM_id, RS_id in product(OMIM_ids, RS_ids):
            info_dict['OMIM_id'] = OMIM_id
            info_dict['RS'] = RS_id
            yield info_dict
    else:
        yield info_dict
    """
    """
    # original:
    for OMIM_id, RS_id in product(OMIM_ids, RS_ids):
        info_dict['OMIM_id'] = OMIM_id
        info_dict['RS'] = RS_id
        yield info_dict
    """

class ClfDatalines:
    """
    HK, 2020-11-21
    """
    def __init__(self, base_string_encoding):
        if base_string_encoding == 'integer':
            self.nuc_encoder = NucEncoder()
            self.base_string_encoder = self._integer_encoding
            self.get_dataline = self._get_dataline_integer_encoding

        elif base_string_encoding == 'base_string_length':
            self.get_dataline = self._get_datalines_base_string_length

    def _integer_encoding(self, base_string: str):
        """
        Encodes nucleobases into integers. Each nucleobase will be encoded into an integer individually.
        Example:
            GATTACA will be encoded to 2033010
        """
        return self.nuc_encoder.encode(None) if str(base_string) == 'nan' else self.nuc_encoder.encode(base_string)

    def _get_dataline_integer_encoding(self, row):
        if (str(row.ref) != 'nan') and (str(row.alt) != 'nan'):
            line = {
                'pos': row['pos'],
                'ref': self.base_string_encoder(row['ref']),
                'alt': self.base_string_encoder(row['alt']),
                'label': 1 if row['CLNSIG'] == 'Pathogenic' else 0
            }
            yield line

    def _get_datalines_base_string_length(self, row):
        yield {'pos': row['pos'], 'label': 1 if row['CLNSIG'] == 'Pathogenic' else 0,
               'length_var': len(row.alt) if str(row.alt) != 'nan' else 0}

    def __call__(self, df_clinvar: pd.DataFrame):
        for idx, row in tqdm(df_clinvar.iterrows(), total=len(df_clinvar), postfix='iterating df_clinvar'):
            yield from self.get_dataline(row)


def df_clinvar_to_clf_data(df_clinvar: pd.DataFrame, base_string_encoding: str = 'integer') -> pd.DataFrame:
    """
    # HK, 2020-11-21
    """
    clf_datalines = ClfDatalines(base_string_encoding=base_string_encoding)
    return pd.DataFrame(data=clf_datalines(df_clinvar))


def clinvar_datalines(vcf: idiva.io.ReadVCF):
    """
    HK, 2020-11-12
    """
    for idx, line in tqdm(enumerate(vcf.datalines), postfix='reading clinvar file'):
        for info_dict in get_info_dict(line.info):
            if info_dict['CLNVC'] == 'single_nucleotide_variant':
                line_dict = {k: line.__dict__[k] for k in line.__dict__.keys() if k != 'info'}
                line_dict = dict(line_dict, **info_dict)
                yield line_dict


def clinvar_to_df(vcf: idiva.io.ReadVCF) -> pd.DataFrame:
    """
    Creates a dataframe from the clinvar file. Adds all the INFO fields as additional columns.
    HK, 2020-11-12
    """

    return pd.DataFrame(data=clinvar_datalines(vcf)).astype({'ref': str, 'alt': str})


def clinvar_rsid_clnsig(which='vcf_37') -> pd.DataFrame:
    """
    RA, 2020-12-05
    """

    def maker():
        from idiva.io.vcf import ReadVCF
        with clinvar_open(which) as fd:
            df = pd.DataFrame(data=(clinvar_datalines(ReadVCF(fd))))
            df = df[["RS", "CLNSIG"]].rename(columns={'RS': "ID", 'CLNSIG': "ClnSig"})
            df = df[df.ID.fillna('').str.contains(r"^rs[0-9]+$")]
            df = df.groupby('ID', as_index=False)
            df = df.agg({'ClnSig': lambda s: F'"{", ".join(sorted(set(map(str, s))))}"'})
            return df

    from idiva.io import cache_df
    return cache_df("rsid_clnsig", [which, 'v1'], df_maker=maker)


if __name__ == '__main__':
    df = clinvar_rsid_clnsig()
    print(df)
