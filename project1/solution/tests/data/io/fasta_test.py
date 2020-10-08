fasta=Fasta("tests\data\io\genome.chr22.5K.fa")
with open("genome_correct.txt","w",encoding="utf-8") as file:
    for i in fasta.yield_sequences(4):
        file.write(str(i)+'\n')