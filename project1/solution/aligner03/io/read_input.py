import io
import re
from collections import namedtuple  


fastq = namedtuple('fastq',['identifier','bases','probs']) 
fasta = namedtuple('fasta', ['bases', 'position'])

def read_fastq_sample(file):
#reads one fq-sample(4 lines) and returns a (identifier,bases,probabilities) tuple or None, if the file is empty, "yielding"
    identifier=file.readline().rstrip()
    if identifier: #check if next line exists
        bases=file.readline().rstrip()
        junk=file.readline() #just the +
        probs=file.readline().rstrip()
        return(fastq(identifier,bases,probs))
    else:
        return None
        
def read_fastq_file(path):
#reads one whole fq-file and returns an array of (identifier,bases,probabilities), slower than read_fastq_sample 
    return_object=[]
    with open(path,"r",encoding="utf-8") as file:
        while True:
            current=read_fastq_sample(file)
            if current:
                #Hendrick put your function here
                return_object.append(current)
            else:
                break
    return return_object
        
def create_first_fasta_string(file,k):
# returns first current string of length k, a buffer from file (may be empty) and if where the first non-N is in the usefull variable
    current_string=''
    while len(current_string)<k:
        buffer=file.readline().rstrip()
        if buffer:
            current_string=current_string+buffer
        else:
            print('The sequence is shorter than k, please adjust k.')
            exit()
    buffer=current_string[k:]
    current_string=current_string[:k]
    
    #check where first not N base is
    match=re.search(r"[ACGTU]",current_string[::-1])
    if match:
        usefull=match.start()
    else:
        usefull=k
    return current_string, list(buffer), usefull
    
def create_buffer(file):
#fetch new buffer from file, return None if file ended/no buffer could be created
    buffer=file.readline().rstrip()
    if not buffer:
        return None
    return list(buffer)
    
def yield_fasta_sequences(current_string, buffer, usefull, k, position):
#yields sequences from current_string and the buffer, might be [] if only N's, also manipulates the usefull variable and gives the correct position
    return_object=[]
    buffer_length=len(buffer)
    while True:
        if usefull<k:
            return_object.append(fasta(current_string, position))
        #book keeping    
        position+=1
        new_char=buffer.pop(0)
        buffer_length-=1
        current_string=current_string[1:]+new_char
        #check if string is not all N's
        if new_char=='N':
            usefull+=1
        else:
            usefull=0
        #one last iteration if buffer is empty
        if  buffer_length<=0:
            if usefull<k:
                return_object.append(fasta(current_string, position))
            #book keeping    
            position+=1
            break
    return return_object, current_string, usefull, position
                    
def read_fasta_file(path,k):
#reads whole fa file on path and returns the sequences of length k, only use as template please
    if k<1:
        print('k must be bigger 0 to make sense.')
        exit()
    return_object=[]
    position=0
    with open(path,"r",encoding="utf-8") as file:
        #title
        return_object.append(file.readline().rstrip()[1:])
        current_string, buffer, usefull=create_first_fasta_string(file,k)
        sequence_list=[]
        #iterate file
        while True:
            sequences, current_string, usefull, position=yield_fasta_sequences(current_string, buffer, usefull, k, position)
            #Luca, use sequences object. Might be [] if all N.
            sequence_list.append(sequences)
            
            
            #end if there is no line left
            buffer=create_buffer(file)
            if not buffer:
                break
    return_object.append(sequence_list)
    return return_object
        
if __name__ == '__main__':
# Debug
    print(read_fastq_file('input\data_small\output_tiny_30xCov1.fq'))
    print(read_fasta_file('input\data_small\genome.chr22.5K.fa',5))
    pass