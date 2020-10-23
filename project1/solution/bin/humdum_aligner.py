# RA, 2020-10-13

from argparse import ArgumentParser
from pathlib import Path

from humdum.main import AllTheKingsHorses
from humdum.utils import assert_exists


def get_args():
    parser = ArgumentParser(description="Align reads to a reference genome.")

    parser.add_argument("fasta", type=str, nargs=1, help="Reference genome as FASTA file.")
    parser.add_argument("fastq", type=str, nargs=2, help="Two FASTQ files.")

    args = parser.parse_args()

    return {
        'fa': assert_exists(Path(args.fasta[0])),
        'fq1': assert_exists(Path(args.fastq[0])),
        'fq2': assert_exists(Path(args.fastq[1])),
    }


def main():
    sam = AllTheKingsHorses.from_files(**get_args())

    for header in sam.headers:
        print(header)

    for alignment in sam.alignments:
        print(alignment)


if __name__ == '__main__':
    main()
