# RA, 2020-11-05

import typing

from unittest import TestCase
from pathlib import Path

from idiva.utils import unlist1, at_most_n, first

vcf_file = unlist1((Path(__file__).parent / "data_for_tests/mini_vcf").glob("*.vcf"))


class TestPlainVCF(TestCase):
    def test_import(self):
        from idiva.io.vcf import ReadVCF

    def test_noimport(self):
        with self.assertRaises(ImportError):
            from idiva.io import PlainVCF

    def test_sanity0(self):
        from idiva.io.vcf import ReadVCF
        with open(vcf_file, mode='r') as fd:
            list(ReadVCF(fd))

    def test_sanity1(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf:
            list(vcf)

    def test_sanity2(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf:
            with ReadVCF.open(vcf) as vcf:
                list(vcf)

    def test_rewinds(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf:
            assert isinstance(vcf, ReadVCF)

            with vcf.rewind_when_done:
                a = list(map(str, vcf))

            with vcf.rewind_when_done:
                b = list(map(str, vcf))

            self.assertListEqual(a, b)

    def test_multi_open1(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf_top:
            assert isinstance(vcf_top, ReadVCF)

            with ReadVCF.open(vcf_top) as vcf:
                a = list(map(str, vcf))

            with ReadVCF.open(vcf_top) as vcf:
                b = list(map(str, vcf))

        self.assertListEqual(a, b)

    def test_multi_open2(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf:
            assert isinstance(vcf, ReadVCF)

            a = list(map(str, vcf))

            with ReadVCF.open(vcf) as vcf:
                b = list(map(str, vcf))

        self.assertListEqual(a, b)

    def test_meta_accuracy(self):
        from idiva.io.vcf import ReadVCF
        with open(vcf_file, mode='r') as fd:
            reference = {
                "INFO": {
                    "NS": {"Number": 1, "Type": "Integer", "Description": '"Number of Samples With Data"'},
                    "DP": {"Number": 1, "Type": "Integer", "Description": '"Total Depth"'},
                    "AF": {"Number": None, "Type": "Float", "Description": '"Allele Frequency"'},
                    "AA": {"Number": 1, "Type": "String", "Description": '"Ancestral Allele"'},
                    "DB": {"Number": 0, "Type": "Flag", "Description": '"dbSNP membership, build 129"'},
                    "H2": {"Number": 0, "Type": "Flag", "Description": '"HapMap2 membership"'},
                },
                "FILTER": {
                    "q10": {"Description": '"Quality below 10"'},
                    "s50": {"Description": '"Less than 50% of samples have data"'},
                },
                "FORMAT": {
                    "GT": {"Number": 1, "Type": "String", "Description": '"Genotype"'},
                    "GQ": {"Number": 1, "Type": "Integer", "Description": '"Genotype Quality"'},
                    "DP": {"Number": 1, "Type": "Integer", "Description": '"Read Depth"'},
                    "HQ": {"Number": 2, "Type": "Integer", "Description": '"Haplotype Quality"'},
                },
                "fileformat": "VCFv4.0",
                "fileDate": "20090805",
                "source": "myImputationProgramV3.1",
                "reference": "1000GenomesPilot-NCBI36",
                "phasing": "partial",
            }

            self.assertDictEqual(ReadVCF(fd).meta, reference)

    def test_datalines_accuracy(self):
        from idiva.io.vcf import ReadVCF
        with ReadVCF.open(vcf_file) as vcf:
            reference = [
                "20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.",
                "20	17330	None	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3",
                "20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4",
                "20	1230237	None	T	None	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2",
                "20	1234567	microsat1	GTCT	G,GTACT	50	PASS	NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3",
            ]
            candidate = list(map(str, vcf))
            self.assertListEqual(reference, candidate)

    def test_dataline_types(self):
        from idiva.io.vcf import ReadVCF, RawDataline
        with ReadVCF.open(vcf_file) as vcf:
            candidate = first(vcf)
            self.assertIsInstance(candidate, RawDataline)

            self.assertIsInstance(candidate.pos, int)
            self.assertIsInstance(candidate.qual, float)
            self.assertIsInstance(candidate.info, str)
