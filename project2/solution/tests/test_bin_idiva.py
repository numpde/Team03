# RA, 2020-10-25
# RA, 2020-12-01

"""
Simulates running the idiva as a script.
"""

from idiva import log

import os
import sys

from pathlib import Path
from unittest import TestCase

from idiva.utils import relpath
from idiva.utils.testing import import_from_source, whatsmyname, redirect

# Set working dir (simulate user)
basepath = Path(__file__).parent.parent
assert str(basepath).endswith("solution")
os.chdir(basepath)

workspace = Path(__file__).parent / "data_for_tests/test_bin/idiva/output/"
workspace.mkdir(exist_ok=True, parents=True)


class TestBin(TestCase):
    def test_idiva_small(self):
        mod = import_from_source("main", str(basepath / "bin/main.py"))

        (case, ctrl) = sorted(map(relpath, Path.cwd().glob("../input/head/*_v2.vcf")))

        out_dir = workspace / whatsmyname()
        out_dir.mkdir(exist_ok=True)

        sys.argv = ["main.py", str(relpath(case)), str(relpath(ctrl)), str(out_dir)]

        with redirect(out_dir / "stdout.log"):
            assert hasattr(mod, 'main')
            mod.main()

    def test_idiva_large(self):
        mod = import_from_source("main", str(basepath / "bin/main.py"))

        (case, ctrl) = sorted(map(relpath, Path.cwd().glob("../input/full/*_v2.vcf.gz")))

        out_dir = workspace / whatsmyname()
        out_dir.mkdir(exist_ok=True)

        sys.argv = ["main.py", str(relpath(case)), str(relpath(ctrl)), str(out_dir)]

        with redirect(out_dir / "stdout.log"):
            assert hasattr(mod, 'main')
            mod.main()
