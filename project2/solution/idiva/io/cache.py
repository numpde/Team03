# RA, 2020-11-22

from idiva import log

import pathlib
import pandas
import hashlib
import base64
import typing

BASE = pathlib.Path(__file__).parent.parent / "_cache"
BASE.mkdir(parents=False, exist_ok=True)

SEP = "\t"


def cache_df(
        name: str,
        key: typing.Union[str, typing.List[str]],
        df_maker: typing.Callable[[], pandas.DataFrame],
) -> pandas.DataFrame:
    """
    Usage:
        df = cache_df(name="my_dataframe", key="v0.1", df_maker=fn_that_makes_df)

    `name` should be a human-readable mnemonic that will appear as part of the file name.
    `key` is a string or list-of-strings that will be hashed and appended to file name.
    `df_maker` is a function that takes no arguments and returns a pandas DataFrame.

    RA, 2020-11-22
    """

    assert name
    assert key is not None, "Provide an empty list or an empty string as `key`."
    assert df_maker is not None

    if isinstance(key, list):
        key = "+".join(key)

    assert isinstance(key, str)

    key = base64.urlsafe_b64encode(hashlib.sha256(key.encode()).digest()).decode()

    name = name + (("__" + key[0:12]) if key else "")
    file = (BASE / name).with_suffix(".gz")
    assert not file.is_dir()

    if True:
        # Fix legacy filename
        old_name = name + (("_" + key[0:12]) if key else "")
        old_name = (BASE / old_name).with_suffix(".gz")
        if old_name.is_file():
            import os
            os.rename(old_name, file)

    if file.is_file():
        log.debug(F"Loading DataFrame from {file.name}.")
        df = pandas.read_csv(file, sep=SEP, compression="infer")
        df.set_index(df.columns[0], inplace=True)
        # https://stackoverflow.com/a/61390345
    else:
        log.info(F"Creating a potentially large DataFrame for the first time ({file.name}).")
        df = df_maker()
        log.debug(F"Saving DataFrame to {file.name}.")
        df.to_csv(file, sep=SEP, compression="gzip", index=True)

    return df
