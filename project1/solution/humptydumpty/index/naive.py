# RA, 2020-10-10


class NaiveIndex:
    def __init__(self, ref):
        self.ref = ref

    def query(self, kmer):
        from aligner03.utils import find_all
        return find_all(pattern=kmer, template=self.ref)

    def __len__(self):
        return len(self.ref)

    def __str__(self):
        return str(self.ref)
