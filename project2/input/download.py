# RA, 2020-10-05

import io
from pathlib import Path
from tcga.utils import download

URLS = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control.vcf",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed.vcf",
}

CACHE = Path(__file__).parent / "download_cache"
download = download.to(abs_path=CACHE)

HEAD = Path(__file__).parent / "head"
HEAD.mkdir(parents=True, exist_ok=True)

# Number of datalines for the `head` preview of VCF
N = 100

for url in URLS.values():
    data = download(url).now

for k in URLS:
    data = download(URLS[k]).now
    vcf: io.FileIO
    with data.open(mode='r') as vcf:
        head = HEAD / Path(data.meta['source']).name
        assert (head.suffix == ".vcf")
        with head.open(mode='w') as fd:
            # Write meta and header
            line = "##"
            while line.startswith("##"):
                line = vcf.readline().strip()
                print(line, file=fd)
            # Write a few datalines
            for _ in range(N):
                line = vcf.readline().strip()
                print(line, file=fd)
