from pathlib import Path
from unittest import TestCase

from idiva.db.dbSNP import dbSNP_to_df
from idiva.io.vcf import ReadVCF

mini_vcf = (Path(__file__).parent) / "data_for_tests/mini_dbSNP/mini_GRCh37_latest_dbSNP_all.vcf"


class Test_dbSNP(TestCase):
    def test_create_df(self):

        with open(mini_vcf, mode='r') as fd:
            df = dbSNP_to_df(ReadVCF(fd))

        assert len(df)

