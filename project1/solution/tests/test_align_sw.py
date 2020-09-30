from align import Smith_waterman_aligner


def test_smith_waterman_aligner():
    test_waterman = Smith_waterman_aligner(ref='ATGGCCTC', query='ACGGCTC', gap_cost=4, match_score=1, mismatch_cost=3)
    start, end = test_waterman.get_start_end()
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
    assert (test_waterman.scoring_matrix == true_scoring_matrix).all
