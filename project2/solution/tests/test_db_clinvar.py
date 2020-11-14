# RA, 2020-11-11

import io
import typing

from unittest import TestCase
from idiva.utils import at_most_n


class TestClinVar(TestCase):
    def test_howto(self):
        from idiva.db import clinvar_open
        from idiva.io.vcf import ReadVCF
        with clinvar_open() as fd:
            vcf = ReadVCF(fd)
            vcf.meta
            for dataline in vcf:
                dataline

    def test_import(self):
        from idiva.db import clinvar_open

    def test_clinvar_meta(self):
        from idiva.db import clinvar_meta
        meta = clinvar_meta()
        assert ('source' in meta)
        assert ('datetime' in meta)

    def test_open(self):
        from idiva.db import clinvar_open
        with clinvar_open() as fd:
            self.assertIsInstance(fd, io.TextIOBase)

    def test_open_read_manual(self):
        from idiva.db import clinvar_open
        with clinvar_open() as fd:
            self.assertIsInstance(fd, io.TextIOBase)
            reference = ['##fileformat=VCFv4.1', '##fileDate=2020-11-07', '##source=ClinVar']
            candidate = [fd.readline().strip() for __ in range(3)]
            self.assertListEqual(reference, candidate)

    def test_open_read_vcf_meta(self):
        from idiva.db import clinvar_open
        from idiva.io import ReadVCF
        with clinvar_open(which='vcf_37') as fd:
            vcf = ReadVCF(fd)
            assert not hasattr(vcf, "sample_ids")
            print(vcf.header)
        raise NotImplementedError

    def test_open_read_vcf_datalines(self):
        from idiva.db import clinvar_open
        from idiva.io import ReadVCF
        with clinvar_open(which='vcf_37') as fd:
            vcf = ReadVCF(fd)

            reference = [
                "1	865568	846933	G	A	.	.	ALLELEID=824438;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.10:g.865568G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1",
                "1	865583	972363	C	T	.	.	ALLELEID=959431;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.10:g.865583C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1",
                "1	865628	789256	G	A	.	.	AF_ESP=0.00347;AF_EXAC=0.00622;AF_TGP=0.00280;ALLELEID=707587;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.10:g.865628G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=SAMD11:148398;MC=SO:0001583|missense_variant;ORIGIN=1;RS=41285790",
            ]

            from idiva.io.vcf import RawDataline
            datalines: typing.List[RawDataline]
            datalines = list(at_most_n(vcf, n=len(reference)))

            self.assertIsInstance(datalines[0], RawDataline)

            candidate = list(map(str, datalines))
            self.assertListEqual(reference, candidate)

            self.assertEqual(datalines[0].ref, 'G')
            self.assertEqual(datalines[1].ref, 'C')
            self.assertEqual(datalines[2].ref, 'G')
