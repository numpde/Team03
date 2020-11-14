# RA, 2020-11-14

"""
Checks certain assumptions about the VCF files.
"""

import typing
from idiva.io.vcf import ReadVCF, RawDataline, is_genomic_string


class bag_of_assumptions:
    @classmethod
    def samples_column(cls, fd):
        from idiva.io.vcf import parse_gt
        for dataline in ReadVCF(fd):
            for gt in dataline.samples:
                try:
                    (a, b) = parse_gt(gt)
                except:
                    raise RuntimeError(F"Could not parse genotype: {gt}")

    @classmethod
    def ref_alt_columns(cls, fd):
        vcf = ReadVCF(fd)
        special = {F"<{k}>" for k in vcf.meta['ALT'].keys()}

        for dataline in vcf:
            ref = dataline.ref
            alt = dataline.alt.split(',')
            if is_genomic_string(ref):
                # Cannot assume:
                # assert all(is_genomic_string(a) for a in alt)
                pass
            else:
                assert ref in special
                assert all((a in special) for a in alt)

    @classmethod
    def ref_column(cls, fd):
        vcf = ReadVCF(fd)

        TCGA = {"T", "C", "G", "A"}
        special = {F"<{k}>" for k in vcf.meta['ALT'].keys()}

        for dataline in vcf:
            assert "," not in dataline.ref

            checks = {
                'single nt': dataline.ref in TCGA,
                'multi nt': set(dataline.ref).issubset(TCGA),
                'special': (dataline.ref in special),
            }

            if not any(checks.values()):
                print(F"REF = '{dataline.ref}' does not fit any known format.")
                print(F"ALT = '{dataline.alt}'.")
                raise RuntimeError("Assumption on REF column failed.")

    @classmethod
    def alt_column(cls, fd):
        vcf = ReadVCF(fd)

        TCGA = {"T", "C", "G", "A"}
        special = {F"<{k}>" for k in vcf.meta['ALT'].keys()}

        for dataline in vcf:
            checks = [
                {
                    'single nt': alt in TCGA,
                    'multi nt': set(alt).issubset(TCGA),
                    'special': alt in special,
                }
                for alt in dataline.alt.split(',')
            ]

            if not any(any(c.values()) for c in checks):
                print(F"ALT = '{dataline.alt}' does not fit any known format.")
                print(F"REF = '{dataline.ref}'.")
                raise RuntimeError("Assumption on ALT column failed.")

    @classmethod
    def format_is_gt(cls, fd):
        for dataline in ReadVCF(fd):
            assert (dataline.format == "GT")


def check_all(fd: typing.TextIO):
    from idiva.utils import seek0

    with seek0(fd):
        bag_of_assumptions.samples_column(fd)

    with seek0(fd):
        bag_of_assumptions.ref_column(fd)

    with seek0(fd):
        bag_of_assumptions.alt_column(fd)

    with seek0(fd):
        bag_of_assumptions.format_is_gt(fd)
