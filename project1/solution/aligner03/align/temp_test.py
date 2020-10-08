from smith_waterman import Smith_Waterman
import pysam
from typing import Iterable

Read = pysam.libcalignedsegment.AlignedSegment

def read_sam(file) -> Iterable[Read]:
    """
    `file` can be file name or file descriptor
    """

    if (type(file) is str):
        with open(file, mode='r') as file:
            yield from read_sam(file)
        return

    # https://pysam.readthedocs.io/en/latest/api.html
    with pysam.AlignmentFile(file) as af:
        for read in af.fetch():
            yield read

def test_read_sam(file):
    with open(file, mode='rb') as fd_sam:
        for read in read_sam(fd_sam):
            read: Read

            # https://en.wikipedia.org/wiki/SAM_(file_format)
            # Col   Field	Type	Brief description
            #   1   QNAME	String  Query template NAME
            #   2   FLAG	Int	    bitwise FLAG
            #   3   RNAME	String  References sequence NAME
            #   4   POS     Int     1- based leftmost mapping POSition
            #   5   MAPQ	Int     MAPping Quality
            #   6   CIGAR	String	CIGAR String
            #   7   RNEXT	String	Ref. name of the mate/next read
            #   8   PNEXT	Int     Position of the mate/next read
            #   9   TLEN	Int     observed Template LENgth
            #  10	SEQ     String  segment SEQuence
            #  11	QUAL	String  ASCII of Phred-scaled base QUALity+33

            print("QNAME ", read.query_name)
            print("FLAG  ", read.flag)
            print("RNAME ", read.reference_id)
            print("POS   ", read.reference_start)
            print("MAPQ  ", read.mapping_quality)
            print("CIGAR ", read.cigartuples)
            print("RNEXT ", read.next_reference_id)
            print("PNEXT ", read.next_reference_start)
            print("TLEN  ", read.template_length)
            print("SEQ   ", read.query_sequence)
            print("QUAL  ", read.query_qualities)

            # print(read.get_aligned_pairs())
            print(read.get_blocks())

def temp_test():
    from pathlib import Path
    from Bio import SeqIO

    fa = Path(__file__).parent.parent.parent.parent / "input/data_small/genome.chr22.5K.fa"

    template = str(SeqIO.read(fa, format='fasta').seq)

    in_file = list((Path(__file__).parent.parent.parent.parent / "input/data_small/").glob("*.sam")).pop()
    with open(in_file, mode='rb') as fd_sam:
        for read in read_sam(fd_sam):
            read: Read
            ref = template
            query = read.query_sequence
            aligner = Smith_Waterman()
            for alignment in aligner(query=query, ref=ref):
                print(alignment.cigar_string, ' vs ', read.cigarstring)
                print(read.query_qualities, ' vs ', alignment.score)
                x, y, z = alignment.visualize(ref, query)
                print(x)
                print(y)
                print(z)
                print(alignment.get_matching_blocks(), ' vs ', read.cigar)
                assert alignment.cigar_string == read.cigarstring, f'{alignment.cigar_string} is not equal to cigar from sam file {read.cigarstring}'




if __name__ == '__main__':
    temp_test()
