# Jannik G, 2020-10

import os
import sys
class Fasta:
#opened fasta file to yield sequences
    file=None
    title="String"
	
    position=0 #only used for book keeping
    buffer=""
    N_length=0 #length of N's at the back
    def __init__(self, path):
   #initializes the Fasta object with opening the file to-be-read
        try:
            self.file=open(path,"r",encoding="utf-8")
        except FileNotFoundError:
            sys.exit("Could not open the following file: "+str(path))
        self.title=self.file.readline().rstrip()
        self.position=0
    def enlarge_buffer(self):
    #reads a line from the file and adds it to the buffer, returns True if file ended
        read=self.file.readline().rstrip()
        if len(read)>0:
            self.buffer+=read
            return False
        else:
            print("The file at the following path has ended: "+str(self.file.name))
            return True
    def yield_sequences(self, sequence_length):
    #yields (not all N sequence, 0-indented position) or returns None if the file ended
        if not isinstance(sequence_length,int) or sequence_length<1:
            self.file.close()
            sys.exit("Sequence length must be positive integer.")        
        while True:    
            while len(self.buffer)<sequence_length:
                if self.enlarge_buffer():
                    self.file.close()
                    return None            
            #first N_length
            if self.position==0:
                                                #last index not N
                self.N_length=sequence_length-max(self.buffer.rfind("A"),self.buffer.rfind("C"),self.buffer.rfind("T"),self.buffer.rfind("G"),self.buffer.rfind("U"))
            if self.buffer[sequence_length-1]=="N":
                self.N_length+=1
            else:
                self.N_length=0
            #don't yield if only Ns
            if self.N_length<sequence_length:
                yield (self.buffer[:sequence_length],self.position)
            self.buffer=self.buffer[1:]
            self.position+=1