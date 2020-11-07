# RA, 2020-11-05

from unittest import TestCase
from pathlib import Path

from idiva.utils import unlist1, at_most_n, first

vcf_file = unlist1((Path(__file__).parent / "data_for_tests/mini_vcf").glob("*.vcf"))


class TestReaderOfVCF(TestCase):
    def test_import(self):
        from idiva.io.vcf import proxy

    def test_noimport(self):
        with self.assertRaises(ImportError):
            from idiva.io import proxy

    def test_call(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            vcf_proxy(fd)

    def test_reads_header(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            proxy = vcf_proxy(fd)
            header = "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003"
            self.assertEqual(header, proxy.header)

    def test_reads_meta(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            proxy = vcf_proxy(fd)
            reference = ["fileformat=VCFv4.0", 'FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">']
            candidate = [proxy.meta[0], proxy.meta[-1]]
            self.assertListEqual(reference, candidate)

    def test_needs_open_file(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            proxy = vcf_proxy(fd)
        with self.assertRaises(ValueError):
            list(next(iter(proxy.datalines)))

    def test_reads_datalines(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            proxy = vcf_proxy(fd)
            reference = [
                "20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.",
                "20	17330	None	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3",
                "20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4",
            ]
            candidate = list(map(str, at_most_n(proxy.datalines, len(reference))))
            self.assertListEqual(reference, candidate)

    def test_reads_dataline_number(self):
        from idiva.io.vcf import proxy as vcf_proxy
        with open(vcf_file, mode='r') as fd:
            proxy = vcf_proxy(fd)
            self.assertEqual(5, len(list(proxy.datalines)))
