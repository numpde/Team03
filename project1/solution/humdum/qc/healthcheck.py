# RA, 2020-10-25

import pathlib


def healthcheck_index(Index):
    genome_file = pathlib.Path(__file__).parent / "healthcheck_index.fa"
    index = Index.read_or_make(path_to_genome=genome_file)
    kmer = "ACCNAANTCGGCNTG"
    loc = [68]
    assert (index.query(kmer) == loc), "Query failed."
