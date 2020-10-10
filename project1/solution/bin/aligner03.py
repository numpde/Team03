
# TODO:
# Take input filenames as arguments and output SAM file to stdout
# python humdum path/to/genome path/to/1.fq path/to/2.fq > alignment.sam

import sys
import time
import resource

#get input args
input_args=sys.argv
if len(input_args)< 4:
    print('''Please give the arguments in the form:
python humdum.py alignment.sam path/to/genome [path/to/1.fq]+ ''')
    exit()
PATH=input_args.pop(0)
OUTPUT_PATH=input_args.pop(0)
GENOME_PATH=input_args.pop(0)
FASTQ_PATHS=[]

for fastq in input_args:
    FASTQ_PATHS.append(fastq)

#construct genome
t=time.time()
read_fasta(GENOME_PATH)
print("FASTA_time: ", (time.time() - t))

#check memory usage
print("FASTA_mem:", (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


for file in FASTQ_PATHS:
    #construct genome
    t=time.time()
    read_fasta(time)
    print("FASTQ_time:", (time.time() - t))
    
    #check memory usage
    print("FASTQ_mem: ", (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))