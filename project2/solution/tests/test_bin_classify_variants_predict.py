from unittest import TestCase
from idiva.bin.classify_variants import predict
from pathlib import Path


class TestClf(TestCase):
    def test_predict_does_something(self):
        vcf_file_path = Path("tests/data_for_tests/large_head/case_processed.vcf")
        which_clf = 'dummy'
        results = predict(vcf_file_path, which_clf)
        self.assertTrue(len(results))
