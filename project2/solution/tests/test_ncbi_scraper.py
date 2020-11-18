from unittest import TestCase
from idiva.db.ncbi_scraper import get_ncbi_info
from pathlib import Path
from idiva.io import ReadVCF
from idiva.utils import at_most_n

BASE = (Path(__file__).parent) / "data_for_tests/large_head"

PATHS = {
    'ctrl': BASE / "control_v2.vcf",
    'case': BASE / "case_processed_v2.vcf",
}


class NcbiScrapeTest(TestCase):
    def test_get_info(self):
        with open(PATHS['ctrl'], mode='r') as fd:
            for elem in at_most_n(ReadVCF(fd), 5):
                info = get_ncbi_info(elem.id)

        self.assertTrue(info)
