def map_pair(self, read1, read2) -> Iterable[AlignedSegment]:
    (read1, options1) = self.map_one(read1).popitem()
    (read2, options2) = self.map_one(read2).popitem()

    if (read1.is_forward == read2.is_forward):
        raise UnmappedReadpair({'reads': [read1, read2], 'reason': "Directionality"})

    # ...

    read2seg = {}

    for (read, options) in zip([read1, read2], [options1, options2]):
        (i, _, _, j) = self.select_option(options)
        (a, b) = propose_window(read_length=len(read), read_loc=i, ref_loc=j, ...)
        genome_segment = self.ref_genome.seq[a:b]
        alignment = first(self.align(ref=genome_segment, query=read.seq))  # ALIGN!
        loc_in_ref = (alignment.loc_in_ref + a)
        # ... copy from `alignment` to SAM-friendly `seg` and then:
        read2seg[read] = seg

    # Get position of mate
    read2seg[read1].pnext = read2seg[read2].pos
    read2seg[read2].pnext = read2seg[read1].pos

    # ...

    # Yield in the correct order
    for read in [read1, read2]:
        yield read2seg[read]
