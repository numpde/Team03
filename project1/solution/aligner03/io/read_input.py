import io
import re

#Luca dummy implementation for BWT
def BWT(current_string,k, current_position):
    return current_string

def read_fastq_file(path):
    return_object=[]
    with open(path,"r",encoding="utf-8") as file:
        while True:
            current=read_fastq_sample(file)
            if current:
                return_object.append(current)
            else:
                break
    return return_object

def read_fastq_sample(file):
    identifier=file.readline().rstrip()
    if identifier: #check if next line exists
        bases=file.readline().rstrip()
        junk=file.readline() #just the +
        probs=file.readline().rstrip()
        return((identifier,bases,probs))
    else:
        return None
    
def read_fasta(path,k):
    return_object=[]
    current_position=0
    current_string=''
    usefull=0
    with open(path,"r",encoding="utf-8") as file:
        #title
        return_object.append(file.readline().rstrip()[1:])
        
        #create first string and buffer
        while len(current_string)<k:
            buffer=file.readline().rstrip()
            if buffer:
                current_string=current_string+buffer
            else:
                print('The sequence is shorter than k, please adjust.')
                exit()
        buffer=current_string[k:]
        current_string=current_string[:k]
        
        #check where first not N base is
        match=re.search(r"[ACGTU]",current_string[::-1])
        if match:
            usefull=match.start()
        else:
            usefull=k
        #iterate file
        while True:
            #go through buffer
            buffer=list(buffer)
            buffer_length=len(buffer)
            while buffer_length>0:
                if usefull<k:
                    print(BWT(current_string,k, current_position)) #Luca
                #book keeping    
                current_position+=1
                new_char=buffer.pop(0)
                buffer_length-=1
                current_string=current_string[1:]+new_char
                #check if string is not all N's
                if new_char=='N':
                    usefull+=1
                else:
                    usefull=0
            #end if there is no line left
            buffer=file.readline().rstrip()
            if not buffer:
                break
    #clean-up last base
    if usefull<k:
        print(BWT(current_string,k, current_position)) #Luca
    print(return_object)
        
        
if __name__ == '__main__':
    #print(read_fastq_file('input\data_small\output_tiny_30xCov1.fq'))
    #read_fasta('input\data_small\genome.chr22.5K.fa',5)
    pass