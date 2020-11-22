# RA, 2020-11-22

import pandas as pd
from pathlib import Path
from unittest import TestCase


def maker() -> pd.DataFrame:
    df = pd.DataFrame(
        data={'Latin': ["A", "B"], 'Greek': ["Alpha", "Beta"]},
        index=[1, 2],
    )
    return df


class no_new_files:
    def __init__(self, path: Path, pattern: str):
        self.path = path
        self.pattern = pattern

    def __enter__(self):
        self.files = set(self.path.glob(self.pattern))

    def __exit__(self, exc_type, exc_val, exc_tb):
        import os
        for f in set(self.path.glob(self.pattern)).difference(self.files):
            os.remove(f)


class TestCacheDf(TestCase):
    def test_import(self):
        from idiva.io import cache_df

    def test_no_new_files(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as base:
            base = Path(base)
            assert base.is_dir()

            file1 = base / "file1"
            file2 = base / "file2"

            with open(file1, mode='w'):
                pass

            assert file1.is_file() and not file2.is_file()

            with no_new_files(base, "file*"):
                with open(file2, mode='w'):
                    pass
                assert file1.is_file() and file2.is_file()

            assert file1.is_file() and not file2.is_file()

        assert not base.is_dir()


    def test_first_make(self):
        key = "test_read_write_01"
        from idiva.io import cache_df
        from idiva.io.cache import BASE
        with no_new_files(BASE, "*.gz"):
            df_reference = maker()
            df_candidate = cache_df(key, maker)

            self.assertTrue(df_reference.equals(df_candidate))

    def test_make_save_read(self):
        key = "test_read_write_02"
        from idiva.io import cache_df
        from idiva.io.cache import BASE
        with no_new_files(BASE, "*.gz"):
            df_reference = maker()

            def phoney_maker():
                raise RuntimeError

            with self.assertRaises(RuntimeError):
                cache_df(key, phoney_maker)

            self.assertTrue(df_reference.equals(cache_df(key, maker)))
            self.assertTrue(df_reference.equals(cache_df(key, phoney_maker)))
            self.assertTrue(df_reference.equals(cache_df(key, phoney_maker)))
