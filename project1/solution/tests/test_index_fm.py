#from 

from ... import FM_Index

from unittest import TestCase, 

class TestIndex(TestCase):
    def test_constructor():
        pass
    
    def test_dollar():
        pass
    
    def test_perfect_match():
        ref = "TAGAGAGATCGATCGACTGACTGACTCAG"
        query = "ACT"

        fm_index = FM_Index(ref_genome)
        hits = fm_index.query(query)

        for i in hits:
            self.assertEqual(ref[i:len(query)], query)
    
        # Make sure that _all_ occurrances are returned
    
        def find(pattern, template):
            """
            Find all occurrences of pattern string in template string.
            Note: returns a list of starting positions with 1- based indexing.
            """
            import re
            return [(m.start() + 1) for m in re.finditer("(?=" + pattern + ")", template)]    
    

if __name__ == "__main__":
    ref_genome = 'TAGAGAGATCGATCGACTGACTGACTCAG$'

    sample = 'ACTACTGTCA'

    print("reference genome: ", ref_genome)
    print("sample: ", sample, "\n")

    k = 3
    index = FM_Index(ref_genome, k)

    print('kmer match  (', sample[:k], ')')
    match = query_index(sample, index)
    print(match, '\n')

    for m in match:
        assert (ref_genome[m:m+k] == sample[:k])

    print('perfect match')
    perfect_match = perfect_match(sample, index)

    print(perfect_match)
