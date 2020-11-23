# RA, 2020-11-23

import pandas

URL = "http://easybioai.com/sc2disease/static/allgwas.txt"


def allgwas() -> pandas.DataFrame:
    from idiva.download import download
    with download(URL).now.open() as fd:
        return pandas.read_csv(fd, sep='\t', header=None)



if __name__ == '__main__':
    print(allgwas())
