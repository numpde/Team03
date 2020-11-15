# RA, 2020-11-11

import contextlib
import gzip
import io
import typing

import pandas as pd

import idiva.utils
from idiva.io.vcf import ReadVCF
from idiva.utils import at_most_n
from tqdm import tqdm

URL = {
    'vcf_37': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
    'vcf_38': "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
}


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


def get_info_dict(info) -> dict:
    info_dict = {}
    for elem in info.split(';'):
        k, v = elem.split('=')
        info_dict[k] = v
    return info_dict


def clinvar_to_df(output_path, which='vcf_37', make_checkpoints=False) -> pd.DataFrame:
    """
    Creates a dataframe from the clinvar file. Adds all the INFO fields as additional columns.
    """
    df = pd.DataFrame()
    with clinvar_open(which=which) as fd:
        vcf = ReadVCF(fd)
        for idx, line in tqdm(enumerate(vcf.datalines), postfix='reading clinvar file'):
            info_dict = get_info_dict(line.info)
            line_dict = line.__dict__
            del line_dict['info']
            line_dict = dict(line_dict, **info_dict)
            df = df.append(line_dict, ignore_index=True)
            if idx % 500 == 0 and make_checkpoints:
                df.to_csv(output_path, compression='gzip', index=False)
    df.to_csv(output_path, compression='gzip', index=False)
    return df


def clinvar_rs_ids(which='vcf_37'):
    from idiva.io.vcf import ReadVCF
    with clinvar_open(which) as fd:
        for dataline in ReadVCF(fd):
            pass


if __name__ == '__main__':
    with clinvar_open('vcf_37') as fd:
        reader = ReadVCF(fd)
        print(*at_most_n(reader.datalines, n=10), sep='\n')
