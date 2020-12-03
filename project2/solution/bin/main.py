# RA, 2020-12-01

from idiva import log

import io
import argparse

from pathlib import Path
from idiva.io import ReadVCF


def main():
    process(**parse_args())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('case_vcf', type=Path, help="VCF file (cases).")
    parser.add_argument('ctrl_vcf', type=Path, help="VCF file (controls).")
    parser.add_argument('out_dir', type=Path, help="Output folder.")
    parser = parser.parse_args()

    assert isinstance(parser.case_vcf, Path)
    assert isinstance(parser.ctrl_vcf, Path)
    assert isinstance(parser.out_dir, Path)

    assert parser.case_vcf.is_file()
    assert parser.ctrl_vcf.is_file()
    assert not parser.out_dir.is_file()

    return dict(case_vcf=parser.case_vcf, ctrl_vcf=parser.ctrl_vcf, out_dir=parser.out_dir)


def process(*, case_vcf: Path, ctrl_vcf: Path, out_dir: Path):
    from idiva.io import open_maybe_gz
    from idiva.io import head

    with open_maybe_gz(case_vcf) as case_full, open_maybe_gz(ctrl_vcf) as ctrl_full:
        assert isinstance(case_full, io.TextIOBase)
        assert isinstance(ctrl_full, io.TextIOBase)

        with head(case_full) as case_head, head(ctrl_full) as ctrl_head:
            log.info("Processing VCF (HEAD).")
            process_vcf(case=ReadVCF(case_head), ctrl=ReadVCF(ctrl_head), out=(out_dir / "head"))

        log.info("Processing VCF (FULL).")
        process_vcf(case=ReadVCF(case_full), ctrl=ReadVCF(ctrl_full), out=(out_dir / "full"))


def process_vcf(*, case: ReadVCF, ctrl: ReadVCF, out: Path):
    out.mkdir(exist_ok=True, parents=True)

    from idiva.clf.df import c3_df, join
    case_idx = c3_df(case)

    from idiva.stat.vcf_to_fisher import vcf_to_fisher
    df = vcf_to_fisher(case=case, ctrl=ctrl)

    # TODO: pack into INFO format
    (join(case=case_idx, ctrl=df))


if __name__ == '__main__':
    main()
