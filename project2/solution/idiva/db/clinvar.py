# RA, 2020-11-11

import contextlib
import gzip
import io
import typing
import idiva.utils

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


if __name__ == '__main__':
    from idiva.io.vcf import ReadVCF
    from idiva.utils import at_most_n

    with clinvar_open('vcf_37') as fd:
        reader = ReadVCF(fd)
        print(*at_most_n(reader.datalines, 10), sep='\n')
