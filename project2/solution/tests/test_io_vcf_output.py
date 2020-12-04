# RA, 2020-12-03

import io
import sys
from contextlib import redirect_stdout

from unittest import TestCase
from pathlib import Path
from idiva.io.vcf import SEP

BASE = (Path(__file__).parent) / "data_for_tests"

MY_SPACE = BASE / F"my_space/{Path(__file__).stem}"
MY_SPACE.mkdir(exist_ok=True)

PATHS_LARGE_HEAD = {
    'ctrl': BASE / "large_head/control_v2.vcf",
    'case': BASE / "large_head/case_processed_v2.vcf",
}

PATHS_LARGE_FULL = {
    'ctrl': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/control_v2.vcf.gz",
    'case': "https://public.bmi.inf.ethz.ch/eth_intern/teaching/cbm_2020/cbm_2020_project2/case_processed_v2.vcf.gz",
}


class Test_VCF_Printout(TestCase):
    def test_poc_head(self):
        """
        Proof-of-concept.
        """
        from idiva.io import ReadVCF, open_maybe_gz

        for (k, ref_file) in PATHS_LARGE_HEAD.items():
            can_file = (MY_SPACE / F"{sys._getframe().f_code.co_name}__{k}").with_suffix(".log")

            with open_maybe_gz(ref_file, mode='r') as fd_ref:
                assert isinstance(fd_ref, io.TextIOBase)
                vcf = ReadVCF(fd_ref).preload_all()

            with open(can_file, mode='w') as fd_can:
                with redirect_stdout(fd_can):
                    for (k, v) in vcf.meta.items():
                        if isinstance(v, str):
                            if (str(k).lower() == "filedate"):
                                from idiva.io.out import fileDate
                                v = fileDate
                            if (str(k).lower() == "source"):
                                from idiva.io.out import source
                                v = source
                            print(F"##{k}={v}")
                        elif isinstance(v, dict):
                            for (i, v) in v.items():
                                assert isinstance(v, dict)
                                assert v
                                p = ','.join(F"{k}={v if v is not None else '.'}" for (k, v) in v.items())
                                print(F"##{k}=<ID={i},{p}>")

                    print(F"#{SEP.join(vcf.header)}")

                    for dataline in vcf:
                        print(str(dataline))

            from idiva.io import Oneliner
            with open_maybe_gz(ref_file, mode='r') as fd_ref:
                with open_maybe_gz(can_file, mode='r') as fd_can:
                    assert isinstance(fd_ref, io.TextIOBase)
                    assert isinstance(fd_can, io.TextIOBase)
                    lines_ref = list(Oneliner(fd_ref))
                    lines_can = list(Oneliner(fd_can))
                    for (ref, can) in zip(lines_ref, lines_can):
                        if not (ref.startswith("##fileDate") or ref.startswith("##source")):
                            self.assertEqual(can, ref)
                    self.assertEqual(len(lines_can), len(lines_ref))
