# RA, 2020-10-19

from argparse import ArgumentParser
from pathlib import Path

from humdum.qc.coverage import coverage_pbp
from humdum.qc.tlen import tlen_hist

from humdum.utils import relpath
from humdum.utils import Plox


def with_suffix(suffix):
    return (lambda path: Path(str(path) + suffix))


def report_this_file(file):
    print("Check out", relpath(file))
    return file


def get_args():
    parser = ArgumentParser(description="Inspect SAM file for quality metrics.")

    parser.add_argument("sam_file", type=str, nargs=1, help="SAM file.")
    parser.add_argument("output_path", type=str, nargs=1, help="Output path.")

    parsed_args = parser.parse_args()

    class arguments:
        sam_file = Path(parsed_args.sam_file[0])
        output_path = Path(parsed_args.output_path[0])

    assert arguments.sam_file.is_file()
    assert arguments.output_path.is_dir()

    return arguments


def qc1_coverage(sam_file: Path, output_path: Path):
    assert sam_file.is_file()
    assert output_path.is_dir()

    counts = coverage_pbp(sam_file)

    fig_file = with_suffix(".coverage.png")(output_path / sam_file.name)

    with Plox() as px:
        px.a.plot(counts)
        px.a.set_xlabel("Position in genome")
        px.a.set_ylabel("Mapped coverage (counts)")
        px.f.savefig(report_this_file(fig_file))


def qc2_tlenhist(sam_file: Path, output_path: Path):
    assert sam_file.is_file()
    assert output_path.is_dir()

    hist = tlen_hist(sam_file)

    fig_file = with_suffix(".tlenhist.png")(output_path / sam_file.name)

    with Plox() as px:
        px.a.bar(x=hist.length, height=hist.counts, width=1)
        px.a.set_xlabel("Transcript length")
        px.a.set_ylabel("Counts")
        px.f.savefig(report_this_file(fig_file))


def main():
    args = get_args()

    qc1_coverage(args.sam_file, args.output_path)
    qc2_tlenhist(args.sam_file, args.output_path)


if __name__ == '__main__':
    main()

    # Ideas:
    # https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html
    # How many reads map to more than 10 locations on the genome?
    # How many reads are unmapped due to read length?
    # What is the average mapped length per read?
