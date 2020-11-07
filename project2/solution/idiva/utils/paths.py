# RA, 2020-10-09

import os
import pathlib


def relpath(path):
    return os.path.relpath(str(path), os.getcwd())


def assert_exists(file):
    assert pathlib.Path(file).exists(), F"Oh, file `{file}`, where art thou?!"
    return file
