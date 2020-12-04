# RA, 2020-10-25
# RA, 2020-12-01

"""
Simulates running the idiva as a script.
"""

from idiva import log

import os
import sys
import types

from pathlib import Path
from unittest import TestCase

from idiva.utils import relpath, unlist1
from contextlib import redirect_stdout

#
basepath = Path(__file__).parent.parent
assert str(basepath).endswith("solution")
os.chdir(basepath)

out_path = Path(__file__).parent / "data_for_tests/test_bin/idiva/output/"
out_path.mkdir(exist_ok=True, parents=True)


# https://stackoverflow.com/questions/33469246/how-do-i-write-python-unit-tests-for-scripts-in-my-bin-directory
def import_from_source(name: str, file_path: str) -> types.ModuleType:
    from importlib.machinery import ModuleSpec, SourceFileLoader
    from importlib.util import spec_from_loader, module_from_spec
    loader: SourceFileLoader = SourceFileLoader(name, file_path)
    spec: ModuleSpec = spec_from_loader(loader.name, loader)
    module: types.ModuleType = module_from_spec(spec)
    loader.exec_module(module)
    return module


class TestBin(TestCase):
    def test_idiva_small(self):
        mod: types.ModuleType = import_from_source("main", str(basepath / "bin/main.py"))

        (case, ctrl) = sorted(map(relpath, Path.cwd().glob("../input/head/*_v2.vcf")))
        out_dir = out_path / sys._getframe().f_code.co_name

        stdout = out_dir / "stdout.log"
        stdout.parent.mkdir(exist_ok=True, parents=True)

        sys.argv = ["main.py", str(relpath(case)), str(relpath(ctrl)), str(out_dir)]

        with stdout.open(mode="w") as fd:
            with redirect_stdout(fd):
                mod.main()

    def test_idiva_large(self):
        mod: types.ModuleType = import_from_source("main", str(basepath / "bin/main.py"))

        (case, ctrl) = sorted(map(relpath, Path.cwd().glob("../input/full/*_v2.vcf.gz")))
        out_dir = out_path / sys._getframe().f_code.co_name

        stdout = out_dir / "stdout.log"
        stdout.parent.mkdir(exist_ok=True, parents=True)

        sys.argv = ["main.py", str(relpath(case)), str(relpath(ctrl)), str(out_dir)]

        with stdout.open(mode="w") as fd:
            with redirect_stdout(fd):
                mod.main()
