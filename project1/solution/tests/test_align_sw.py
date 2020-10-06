from aligner03.align import Smith_waterman_aligner


def test_smith_waterman_aligner(verbose=0):
    """
    Test if sw_aligner finds the right scoring matrix and the right matching blocks
    """
    test_waterman = Smith_waterman_aligner(ref='ATGGCCTC', query='ACGGCTC', gap_cost=4, match_score=1, mismatch_cost=3)
    matching_blocks = test_waterman.get_matching_blocks()
    # true values from http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman#
    true_scoring_matrix = [[0, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 1, 0],
                           [0, 0, 0, 1, 1, 0, 0, 0],
                           [0, 0, 0, 1, 2, 0, 0, 0],
                           [0, 0, 1, 0, 0, 3, 0, 1],
                           [0, 0, 1, 0, 0, 1, 0, 1],
                           [0, 0, 0, 0, 0, 0, 2, 0],
                           [0, 0, 1, 0, 0, 1, 0, 3]]
    # true matching blocks calculated on paper
    true_matching_blocks = [[3, 5]]

    if verbose:
        x, y, z = test_waterman.visualize(test_waterman.ref, test_waterman.query,
                                          test_waterman.compress_cigar(test_waterman.cigar))
        print(x)
        print(y)
        print(z)
        print('matching blocks = ', matching_blocks)

    assert (test_waterman.scoring_matrix == true_scoring_matrix).all
    assert matching_blocks == true_matching_blocks

def test_sw_on_data():
    """
    Test if sw_aligner can find position of read in template reference genome
    """
    from pathlib import Path
    from Bio import SeqIO

    fa = Path(__file__).parent.parent.parent / "input/data_small/genome.chr22.5K.fa"

    template = str(SeqIO.read(fa, format='fasta').seq)

    # This is from output_tiny_30xCov1.fq: 22_5K-1168/1
    r1 = "TCTGGGCCTCCCAACCCTGAGTTTTTATAATAGGCCCCAGGCCAGGTGGTTAACAGAGGTCTGGGGCATTGCAGGGGGACAGAGGAGGACATATGTCCCTATTGGCCATTGTAGAGTCCCTTCCA"

    # String operators
    reverse = (lambda s: ''.join({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[c] for c in s[-1::-1]))

    def find(pattern, template):
        """
        Find all occurrences of pattern string in template string.
        Note: returns a list of starting positions with 1- based indexing.
        """
        import re
        return [(m.start() + 1) for m in re.finditer("(?=" + pattern + ")", template)]

    reversed = reverse(r1)
    true_index = find(reversed, template)[0]
    test_waterman_reversed = Smith_waterman_aligner(ref=template, query=reversed, gap_cost=1, match_score=3,
                                                    mismatch_cost=3)
    beg, end = test_waterman_reversed.get_start_end()
    assert beg == true_index


if __name__ == '__main__':
    test_smith_waterman_aligner(1)
    test_sw_on_data()