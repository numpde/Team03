# RA, 2020-11-22

import pathlib
import pandas
import hashlib
import base64
import typing

BASE = pathlib.Path(__file__).parent.parent / "_cache"
BASE.mkdir(parents=False, exist_ok=True)

SEP = "\t"


def cache_df(
        key: typing.Union[str, typing.List[str]],
        df_maker: typing.Callable[[], pandas.DataFrame],
) -> pandas.DataFrame:
    """
    Usage:
        df = cache_df("v0.1", fn_that_makes_df)

    RA, 2020-11-22
    """

    if isinstance(key, list):
        key = "+".join(key)

    assert isinstance(key, str)

    name = base64.urlsafe_b64encode(hashlib.sha256(key.encode()).digest()).decode()
    file = (BASE / name).with_suffix(".gz")
    assert not file.is_dir()

    if file.is_file():
        df = pandas.read_csv(file, sep=SEP, compression="infer", index_col=0)
    else:
        df = df_maker()
        df.to_csv(file, sep=SEP, compression="gzip", index=True)

    return df
