# RA, 2020-11-23

import pandas

URL = "http://easybioai.com/sc2disease/static/allgwas.txt"


def allgwas() -> pandas.DataFrame:
    from idiva.download import download
    with download(URL).now.open() as fd:
        return pandas.read_csv(fd, sep='\t', names=["gene", "rs", "p", "chrom", "pos", "disease"])


if __name__ == '__main__':
    from collections import Counter
    allgwas = allgwas()

    # https://www.snpedia.com/index.php/Rs7903146
    (rs, __) = max(Counter(allgwas.rs).items(), key=(lambda x: x[1]))
    assert rs == "rs7903146"
    print(allgwas[allgwas.rs == rs])

