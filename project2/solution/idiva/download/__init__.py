# RA, 2020-11-11

import pathlib
import tcga.utils

download = tcga.utils.download.to(abs_path=pathlib.Path(__file__).parent.resolve())
# Usage: data = download(url).now
#        print(data.meta)
#        with data.open() as fd: print(fd.readline())
