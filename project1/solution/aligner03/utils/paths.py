# RA, 2020-10-09

import os


def relpath(path):
    return os.path.relpath(str(path), os.getcwd())
