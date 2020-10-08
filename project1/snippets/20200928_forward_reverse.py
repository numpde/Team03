# RA, 2020-09-28
# How are paired reads oriented in the reference genome?

from pathlib import Path
from Bio import SeqIO

fa = Path(__file__).parent.parent / "input/data_small/genome.chr22.5K.fa"

template = str(SeqIO.read(fa, format='fasta').seq)

# This is from output_tiny_30xCov1.fq: 22_5K-1168/1
r1 = "TCTGGGCCTCCCAACCCTGAGTTTTTATAATAGGCCCCAGGCCAGGTGGTTAACAGAGGTCTGGGGCATTGCAGGGGGACAGAGGAGGACATATGTCCCTATTGGCCATTGTAGAGTCCCTTCCA"
# This is from output_tiny_30xCov2.fq: 22_5K-1168/2
# It is also the same as 22_5K-880/1 and 22_5K-712/2
r2 = "CTGGGTGACAGAAGGAGACTCTGTCTAGGAAAAACAAAAATTTTGCAGAAAATATGTTGAGACAGGATCTATGTTGCCCAGGCTGGTCTTGAACTCCTGGGCTCAAGCAATCCTTCTGCTTCCGA"

# String operators
forward = (lambda s: s)
reverse = (lambda s: ''.join({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[c] for c in s[-1::-1]))


# Expect that one of the paired reads will be matched as reverse-complement
# Note: if read the reverse-complement r' of r is matched, it shows up in the SAM file as r' (not r)

# NOTE: See aligner03.utils.strings
def find(pattern, template):
    """
    Find all occurrences of pattern string in template string.
    Note: returns a list of starting positions with 1- based indexing.
    """
    import re
    return [(m.start() + 1) for m in re.finditer("(?=" + pattern + ")", template)]


for (op, op_name) in [(forward, "forward"), (reverse, "reverse-complement")]:
    print(F"Matches for {op_name} string")
    print("    Read 1:", find(op(r1), template))
    print("    Read 2:", find(op(r2), template))

    # print(op(r1))
    # print(op(r2))
