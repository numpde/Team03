# RA, 2020-10-10
    
from humdum.io import from_fasta
from humdum.utils import unlist1


class NaiveIndex:
    def __init__(self, ref):
        self.ref = ref

    def query(self, kmer):
        from humdum.utils import find_all
        return find_all(pattern=kmer, template=self.ref)

    def __len__(self):
        return len(self.ref)

    def __str__(self):
        return str(self.ref)

    @classmethod
    def read_or_make(cls, *, path_to_genome, ignored=None):
        return cls(unlist1(list(from_fasta(path_to_genome))).seq)
