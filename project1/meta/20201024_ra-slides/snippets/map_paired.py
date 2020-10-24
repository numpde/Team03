class AllTheKingsHorses:
    """
    Attempts to undo shotgunning.
    """

    # ...

    def map_paired(self, file1, file2) -> Iterable[AlignedSegment]:
        """
        Yield aligned segments.
        """

        assert_order_consistency(file1, file2)

        self.unmapped_readpairs = 0

        for (read1, read2) in zip(from_fastq(file1), from_fastq(file2)):
            try:
                yield from self.map_pair(read1, read2)
            except UnmappedReadpair as ex:
                self.unmapped_readpairs += 1
