# RA, 2020-10-25

"""
This simulates running the aligner as a script,
and verifies the md5 checksum of the SAM output
for data_small as input.
"""

import os
import sys
import types
import hashlib

from pathlib import Path
from unittest import TestCase

from humdum.utils import relpath, unlist1
from contextlib import redirect_stdout

path = Path(__file__).parent.parent
assert str(path).endswith("solution")
os.chdir(path)


# https://stackoverflow.com/questions/33469246/how-do-i-write-python-unit-tests-for-scripts-in-my-bin-directory
def import_from_source(name: str, file_path: str) -> types.ModuleType:
    from importlib.machinery import ModuleSpec, SourceFileLoader
    from importlib.util import spec_from_loader, module_from_spec
    loader: SourceFileLoader = SourceFileLoader(name, file_path)
    spec: ModuleSpec = spec_from_loader(loader.name, loader)
    module: types.ModuleType = module_from_spec(spec)
    loader.exec_module(module)
    return module


mod: types.ModuleType = import_from_source("humdum_aligner", str(path / "bin/humdum_aligner.py"))


class TestBin(TestCase):
    def test_humdum_aligner_data_small(self):
        fa = unlist1(list(Path.cwd().glob("../input/data_small/*.fa")))
        (fq1, fq2) = sorted(map(relpath, Path.cwd().glob("../input/data_small/*.fq")))

        out = Path.cwd() / "../output/data_small/alignment.sam"
        out.parent.mkdir(exist_ok=True, parents=True)

        sys.argv = ["humdum_aligner.py", str(fa), str(fq1), str(fq2)]

        with out.open(mode="w") as fd:
            with redirect_stdout(fd):
                mod.main()

        # Calculated with md5sum command
        reference = "91832e10cd90a683aad9bc6a58ddda2c"

        with out.open(mode='rb') as fd:
            candidate = hashlib.md5(fd.read()).hexdigest()

        self.assertEqual(reference, candidate)
