# HK, 2020-12-02

from idiva import log

import re
import typing
from collections.abc import MutableMapping
import gzip
import pandas as pd
from pathlib import Path
import idiva.utils
from idiva.io.vcf import ReadVCF

MAX_LEN_DF = 19260539
TOTAL_LEN = 717580549

DTYPES: dict = {'ref': str, 'alt': str, 'chrom': str, 'id': str, 'GNO': bool, 'COMMON': bool, 'INT': bool, 'pos': int,
                'CLNSIG': 'category', 'CLNORIGIN': 'category', 'R5': bool, 'R3': bool, 'DSS': bool, 'ASS': bool,
                'U5': bool,
                'U3': bool,
                'SYN': bool, 'NSN': bool, 'NSM': bool, 'NSF': bool, 'SAO': float, 'FREQ': str, 'VC': str}

dbSNP_URL = "https://www.dropbox.com/s/u2vfggr70prk3wx/GRCh37_latest_dbSNP_all_chrom17.csv.gz?dl=1"


def flatten(d, parent_key='', sep='_') -> dict:
    """
    Flattens a nested dictionary
    Example:
        for {'a':{'b':1, 'c':2}}
        returns {a_b:1, a_c:2}

    HK, 2020-12-02
    """
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def get_info_dict(info: str) -> typing.Iterable[dict]:
    """
    Yields info dict for every RS id.

    HK, 2020-12-02
    """
    RS_ids = [None]
    # freq_dict = {'FREQ': {}}
    info_dict = {}
    # go through info which is semicolon separated
    for elem in info.split(';'):
        try:
            if len(elem.split('=')) == 2:
                # spit into key, value
                k, v = elem.split('=')
                if k in DTYPES:
                    info_dict[k] = int(re.findall(r"[0-9]+", v)[0]) if k in ['CLNSIG', 'CLNORIGIN'] else v

                    # if k == 'FREQ':
                    #     # the FREQ info resembles this: FREQ=KOREAN:0.9891,0.0109|SGDP_PRJ:0,1
                    #     # split first by | then by :
                    #     for freq in v.split('|'):
                    #         country, frequencies = (freq.split(':')[0], freq.split(':')[1])
                    #         freq_dict['FREQ'][country] = frequencies

            elif len(elem.split('=')) == 1 and elem in DTYPES:
                info_dict[elem] = True

        except Exception as e:
            print(e)

    for RS_id in RS_ids:
        info_dict['RS'] = RS_id
        # flat = flatten(freq_dict)
        # info_dict.update(flat)
        yield info_dict


def dbSNP_datalines(vcf: idiva.io.ReadVCF, which_chrom: str = 'NC_000017.10') -> typing.Iterable[dict]:
    """
    HK, 2020-12-02
    """
    from tqdm import tqdm
    from itertools import product
    for idx, line in tqdm(enumerate(vcf.datalines), postfix='reading dbSNP file',
                          total=MAX_LEN_DF):
        if line.chrom == which_chrom:
            for info_dict in get_info_dict(line.info):
                if info_dict['VC'] == 'SNV':
                    # making pos,ref,alt unique:
                    for ref, alt in product(line.ref.split(','), line.alt.split(',')):
                        line_dict = {k: line.__dict__[k] for k in line.__dict__.keys() if k != 'info'}
                        line_dict['alt'] = alt
                        line_dict['ref'] = ref
                        line_dict = dict(line_dict, **info_dict)
                        yield line_dict


def dbSNP_to_df(vcf: idiva.io.ReadVCF) -> pd.DataFrame:
    """
    Creates a dataframe from the clinvar file. Adds all the INFO fields as additional columns.
    HK, 2020-12-02
    """
    # fill nans for type conversion
    try:

        df = pd.DataFrame(data=dbSNP_datalines(vcf))
        dtypes = {k: v for k, v in DTYPES.items() if k in df.columns}
        df = df.fillna(
            value={'GNO': False, 'COMMON': False, 'INT': False, 'R5': False, 'R3': False, 'DSS': False, 'ASS': False,
                   'U5': False, 'U3': False,
                   'SYN': False, 'NSN': False, 'NSM': False, 'NSF': False, 'CLNSIG': 1, 'CLNORIGIN': 0}).astype(dtypes)
    except:
        df = pd.DataFrame(data=dbSNP_datalines(vcf))
    return df


def get_dbSNP_df() -> pd.DataFrame:
    """
    Returns the dbSNP for chrom 17 as dataframe.
    Downloads it if not found, loads it otherwise
    """
    from idiva.download import download
    import gzip

    log.info("Downloading dbSNP excerpt.")
    with download(dbSNP_URL).now.open(mode='rb') as fd:
        with gzip.open(fd, mode='r') as fd:
            df = pd.read_csv(fd)

    return df


def create_dbSNP_df(dbSNP_file_path: str, out_base: Path):
    """
    Converts the dbSNP vcf file to a dataframe
    """
    out_path = out_base / 'GRCh37_latest_dbSNP_all_chrom17.csv.gz'
    print(out_path)
    assert out_base.exists()

    with open(dbSNP_file_path, mode='r') as fd:
        df = dbSNP_to_df(ReadVCF(fd))

    df.to_csv(out_path, index=False, compression="gzip")


if __name__ == '__main__':
    from pathlib import Path
    from tqdm import tqdm

    dbSNP_file_path = '/mnt/data/hendrik/db_SNP/GRCh37_latest_dbSNP_chrom17.vcf'
    out_base = Path(__file__).parent.parent.parent / 'data'
    out_path = out_base / 'GRCh37_latest_dbSNP_all_chrom17.csv.gz'
    print(out_path)
    assert out_base.exists()

    with open(dbSNP_file_path, mode='r') as fd:
        df = dbSNP_to_df(ReadVCF(fd))
    df.to_csv(out_path, index=False, compression="gzip")
