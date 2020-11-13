from unittest import TestCase
from argparse import ArgumentParser
from idiva.main import main


class TestClf(TestCase):
    def test_main(self):
        args = {
            'vcf_file_path': "tests/data_for_tests/large_head/case_processed.vcf",
            'which_clf': 'dummy',
        }
        parser = ArgumentParser()
        flags = parser.parse_args([])
        flags.__dict__.update(args)
        results = main(flags)
        self.assertTrue(len(results))
