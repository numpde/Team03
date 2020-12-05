# RA, 2020-10-05

import io
from contextlib import ExitStack
from pathlib import Path
from tcga.utils import download

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}

CACHE = Path(__file__).parent / "download_cache"
download = download.to(abs_path=CACHE)

HEAD = Path(__file__).parent / "head"
HEAD.mkdir(parents=True, exist_ok=True)

# Number of datalines for the `head` preview of VCF
N = 1000

for url in URLS.values():
    data = download(url).now

for k in URLS:
    data = download(URLS[k]).now
    head = HEAD / Path(data.meta['source']).name

    with ExitStack() as stack:
        src = stack.enter_context(data.open(mode='rb'))

        try:
            import gzip
            src = stack.enter_context(gzip.open(src))
        except:
            raise
        else:
            head = Path(str(head)[:-3])
        finally:
            src = io.TextIOWrapper(src)

        assert (head.suffix == ".vcf")

        with head.open(mode='w') as fd:
            # Write meta and header
            line = "##"
            while line.startswith("##"):
                line = src.readline().strip()
                print(line, file=fd)
            # Write a few datalines
            for _ in range(N):
                line = src.readline().strip()
                print(line, file=fd)
