# RA, 2020-12-01

from idiva import log

import io
import sys
import gzip
import pandas
import argparse

from pathlib import Path
from idiva.io import ReadVCF
from idiva.utils import relpath


def main():
    try:
        process(**parse_args())
    except KeyboardInterrupt:
        import logging
        log.info("Aborted. Bye-bye.")
        logging.shutdown()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('case_vcf', type=Path, help="VCF file (cases).")
    parser.add_argument('ctrl_vcf', type=Path, help="VCF file (controls).")
    parser.add_argument('out_dir', type=Path, help="Output folder.")
    parser = parser.parse_args()

    case_vcf = Path(parser.case_vcf)
    ctrl_vcf = Path(parser.ctrl_vcf)
    out_dir = Path(parser.out_dir)

    assert case_vcf.is_file(), F"{case_vcf} not found."
    assert ctrl_vcf.is_file(), F"{ctrl_vcf} not found."
    assert not out_dir.is_file(), F"{out_dir} is not a folder."

    return dict(case_vcf=case_vcf, ctrl_vcf=ctrl_vcf, out_dir=out_dir)


def process(*, case_vcf: Path, ctrl_vcf: Path, out_dir: Path):
    from idiva.io import open_maybe_gz
    from idiva.io import head

    with open_maybe_gz(case_vcf) as case_full, open_maybe_gz(ctrl_vcf) as ctrl_full:
        assert isinstance(case_full, io.TextIOBase)
        assert isinstance(ctrl_full, io.TextIOBase)

        with head(case_full) as case_head, head(ctrl_full) as ctrl_head:
            log.info("======================")
            log.info("Processing VCF (HEAD).")
            log.info("======================")
            post(process_vcf(case=ReadVCF(case_head), ctrl=ReadVCF(ctrl_head), out=(out_dir / "head")))

        log.info("======================")
        log.info("Processing VCF (FULL).")
        log.info("======================")
        post(process_vcf(case=ReadVCF(case_full), ctrl=ReadVCF(ctrl_full), out=(out_dir / "full")))


def process_vcf(*, case: ReadVCF, ctrl: ReadVCF, out: Path):
    from idiva.clf.df import join
    from idiva.io.out import spit_out_vcf_with_extra_info_no_samples
    from contextlib import redirect_stdout

    out.mkdir(exist_ok=True, parents=True)

    with case.rewind_when_done:
        from idiva.clf.df import c5_df
        info_supp = c5_df(case)

    info_meta = []

    # # #

    from idiva.clf.placeholder import placeholder, failure
    from idiva.stat.vcf_to_fisher import vcf_to_fisher
    from idiva.clf.sc2disease import allgwas

    classifiers = [failure, placeholder, vcf_to_fisher, allgwas]

    for classifier in classifiers:
        with case.rewind_when_done:
            log.info(F"=> Invoking the annotation `{classifier.__name__}`.")
            try:
                response = classifier(case=case, ctrl=ctrl)
                info_meta.append(response.info)

                assert set(response.id_cols).issubset(set(info_supp.columns))
                df = response.df[set(response.id_cols) | set(response.info.keys())]

                info_supp = join(case=info_supp, ctrl=df, how="left", on=list(response.id_cols))

                del response
            except KeyboardInterrupt:
                raise
            except Exception as ex:
                log.exception(F"=> Annotation `{classifier.__name__}` failed ({ex}).")
            else:
                log.info(F"=> Annotation `{classifier.__name__}` OK.")

    # # # #

    vcf_out = (out / "results.vcf.gz")
    log.info(F"=> Writing VCF to: {relpath(vcf_out)} .")

    with gzip.open(vcf_out, mode="wt") as fd:
        with redirect_stdout(fd):
            info_meta = {k: i[k] for i in info_meta for k in i}
            spit_out_vcf_with_extra_info_no_samples(case, info_supp, info_meta)
            sys.stdout.flush()

        fd.flush()

    log.info("=> Done.")

    return vcf_out


def post(vcf_file: Path):
    log.info("=> Entering the postprocessing stage.")

    from idiva.stat.vcf_to_fisher import figure_pvalues
    from idiva.io.vcf import SEP

    with ReadVCF.open(vcf_file) as vcf:
        for px in figure_pvalues(vcf):
            file = vcf_file.parent / px.info['name proposal']
            log.info(F"Saving figure and data to {file}.* .")

            px.f.savefig(file.with_suffix(".png"))

            df: pandas.DataFrame = px.info['df']
            df.to_csv(file.with_suffix(".csv"), sep=SEP)


if __name__ == '__main__':
    main()
