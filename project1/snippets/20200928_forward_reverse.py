# RA, 2020-09-28
# How are paired reads oriented in the reference genome?

fa = "/home/ra/repos/Team03/project1/input/data_small/genome.chr22.5K.fa"

from Bio import SeqIO

ref = SeqIO.read(fa, format='fasta')

# 22_5K-1168/1
r1 = "TCTGGGCCTCCCAACCCTGAGTTTTTATAATAGGCCCCAGGCCAGGTGGTTAACAGAGGTCTGGGGCATTGCAGGGGGACAGAGGAGGACATATGTCCCTATTGGCCATTGTAGAGTCCCTTCCA"
# 22_5K-1168/2
r2 = "CTGGGTGACAGAAGGAGACTCTGTCTAGGAAAAACAAAAATTTTGCAGAAAATATGTTGAGACAGGATCTATGTTGCCCAGGCTGGTCTTGAACTCCTGGGCTCAAGCAATCCTTCTGCTTCCGA"

forward = (lambda s: s)
reverse = (lambda s: ''.join({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[c] for c in s[-1::-1]))

print(ref.seq.find(reverse(r1)))
print(ref.seq.find(forward(r2)))

assert (reverse(r1) in ref.seq)
assert (forward(r2) in ref.seq)
