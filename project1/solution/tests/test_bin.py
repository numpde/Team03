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

basepath = Path(__file__).parent.parent
assert str(basepath).endswith("solution")
os.chdir(basepath)

out_path = Path(__file__).parent / "data_for_tests/test_bin/output/data_small"
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
    def test_humdum_aligner_data_small(self):
        mod: types.ModuleType = import_from_source("humdum_aligner", str(basepath / "bin/humdum_aligner.py"))

        fa = unlist1(list(Path.cwd().glob("../input/data_small/*.fa")))
        (fq1, fq2) = sorted(map(relpath, Path.cwd().glob("../input/data_small/*.fq")))

        out = out_path / "alignment.sam"
        out.parent.mkdir(exist_ok=True, parents=True)

        sys.argv = ["humdum_aligner.py", str(relpath(fa)), str(relpath(fq1)), str(relpath(fq2))]

        with out.open(mode="w") as fd:
            with redirect_stdout(fd):
                mod.main()

        # Calculated with md5sum command
        reference = "91832e10cd90a683aad9bc6a58ddda2c"

        with out.open(mode='rb') as fd:
            candidate = hashlib.md5(fd.read()).hexdigest()

        self.assertEqual(reference, candidate)

    def test_humdum_qc_data_small(self):
        mod: types.ModuleType = import_from_source("humdum_aligner", str(basepath / "bin/humdum_qc.py"))

        sam = unlist1(list(out_path.glob("*.sam")))

        sys.argv = ["humdum_qc.py", str(relpath(sam)), str(relpath(out_path))]

        mod.main()
