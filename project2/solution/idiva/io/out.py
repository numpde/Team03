# RA, 2020-12-03

from idiva import log

import datetime
import idiva.io
import pandas
import typing
import tqdm  # do not delete

fileDate = datetime.datetime.utcnow().strftime("%Y%m%d")
source = '"1000GenomesPhase3Pipeline via CompBiomed@ETHZ -- Team03"'


def spit_out_vcf_with_extra_info_no_samples(
        vcf: idiva.io.ReadVCF,
        info_data: pandas.DataFrame,
        info_meta: dict,
):
    from idiva.io.vcf import SEP

    from copy import deepcopy
    meta = deepcopy(vcf.meta)

    assert 'INFO' in meta

    for (k, v) in info_meta.items():
        assert k not in meta['INFO'].keys()
        meta['INFO'][k] = v

    for (k, v) in meta.items():
        if isinstance(v, str):
            if (str(k).lower() == "filedate"):
                v = fileDate
            if (str(k).lower() == "source"):
                v = source
            print(F"##{k}={v}")
        elif isinstance(v, dict):
            for (i, v) in v.items():
                assert isinstance(v, dict)
                assert v
                p = ','.join(F"{k}={v if v is not None else '.'}" for (k, v) in v.items())
                print(F"##{k}=<ID={i},{p}>")

    print(F"#{SEP.join(vcf.header[0:8])}")

    i = info_data.iterrows()

    from idiva.io.vcf import RawDataline
    from tqdm import tqdm

    for (dataline, (__, info)) in tqdm(zip(vcf, i)):
        assert isinstance(dataline, RawDataline)
        assert (dataline.chrom == info.CHROM)

        # Align on the fly
        while True:
            a = dataline.pos
            b = info.POS
            if (a <= b):
                break
            else:
                (__, info) = next(i)

        dataline.format = None
        dataline.samples = []

        if ((dataline.id, dataline.ref, dataline.alt) == (info.ID, info.REF, info.ALT)):
            for c in info_meta:
                if info_meta[c] is not None:
                    t = dict(info_meta[c]).get('Type')
                    if t is not None:
                        if (t == 'Float'):
                            dataline.info += (F";{c}={info[c]:.3e}")
                            continue
                        if (t == 'Integer'):
                            dataline.info += (F";{c}={info[c]}")
                            continue
                    dataline.info += (F";{c}={info[c]}")

        print(str(dataline), flush=True)
