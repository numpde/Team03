def map_pair(self, read1, read2) -> Iterable[AlignedSegment]:
    (read1, options1) = self.map_one(read1).popitem()
    (read2, options2) = self.map_one(read2).popitem()

    # ...


def map_one(self, read, decide=True) -> Dict[Read, List]:
    """
    The primary purpose of this is to find whether
    the read aligns forward or backward to the reference.
    """

    proposals = {
        r: [
            (loc_in_read, None, qual, loc_in_ref)
            for (loc_in_read, kmer, qual) in random_kmers(r, k=..., maxn=...)
            for loc_in_ref in list(self.index.query(kmer))  # CALL TO INDEX!
        ]
        for r in [read, read.reversed]
    }

    if decide:
        # Keep the read with more matches
        return dict([max(proposals.items(), key=(lambda p: len(p[1])))])
    else:
        return proposals
