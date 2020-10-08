import io
from Bio import SeqIO
Input="output_10xCov1.fq"#all three given files
Output=Input[:-2]+"test"
Output=open(Output, "w")
for record in SeqIO.parse(Input,"fastq"):
    Output.write(str(record.id)+"\n")
    Output.write(str(record.seq)+"\n")
    Output.write(str("".join(map(str,record.letter_annotations["phred_quality"])))+"\n")