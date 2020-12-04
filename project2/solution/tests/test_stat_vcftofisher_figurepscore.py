# RA, 2020-12-04

from unittest import TestCase
from pathlib import Path

from idiva.io import ReadVCF
from idiva.utils import unlist1

basepath = Path(__file__).parent / "data_for_tests/test_bin/idiva/output/test_idiva_large"
assert basepath.is_dir()


class FigurePscoreTest(TestCase):
    def test_import(self):
        from idiva.stat.vcf_to_fisher import figure_pvalues

    def test_figure(self):
        from idiva.stat.vcf_to_fisher import figure_pvalues
        for kind in ["head", "full"]:
            workdir = basepath / kind
            vcf = unlist1(workdir.glob("*.vcf.gz"))
            with ReadVCF.open(vcf) as vcf:
                for px in figure_pvalues(vcf):
                    px.f.savefig((workdir / px.info['name proposal']).with_suffix(".png"))
