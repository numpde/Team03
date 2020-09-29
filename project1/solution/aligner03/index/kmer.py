from collections import defaultdict

#k-mer indexing with a dictionary for a reference genome
#Keys are the k-mer's present in the reference genome
#Values are the corresponding indexes in the reference genome
#the query function will search for a match of the first k letters
#this could be adapted to a shifted matching or a special pattern if desired
class kmer_Index(object):
    
    #init an index for the reference genome using k-mer indexing
    #reference_genome should be a string
    def __init__(self, reference_genome , k):
    
        if(k < 1 or k > len(reference_genome)):
    	    raise IndexError(' k has to be greater than 0 and smaller than the length of the reference genome. Make sure reference genome is not empty')
    	      
        self.k = k
        self.index = defaultdict(list)   
        
        #iterate over all k-mer's and add it to the index
        for i in range (len(reference_genome) - k + 1):
            self.index[reference_genome[i:i+k] ].append(i)
           
                    
    #query index for a given sample
    #this searches ONLY FOR A MATCH OF THE FIRST KMER of the sample
    #sample should be a string
    #return a list of all positions in the reference genome
    def query(self, sample):
    
        if(len(sample) < self.k):
            raise ValueError('sample length needs to be greater than k = '+str(self.k))
    
        kmer = sample[:self.k]
        values = self.index[kmer]

        if(len(values)==0):
            print("kmer not in reference genome")
            return []
            
        return values

            
    #for debugging/testing 
    def print(self):
    
        print("k = " ,self.k)
        print(self.index)
        
#searching for a perfect match of the sample in the reference genome
#sample should be string
#reference_genome should be string
def perfect_match(sample, reference_genome, index):

    k = index.k
    positions = []
    
    for i in index.query(sample):
        if sample[k:] == reference_genome[i+k:i+len(sample)]:
            positions.append(i)
   
    return positions
    
    
#searching for a match of the first k letters of the sample in the reference genome
#sample should be string
#reference_genome should be string
def query_index(sample, reference_genome, index):
   
    return index.query(sample)
        
        
if __name__ == "__main__":
      
    ref_genome = 'GCTACCATCTAGAATCTG'
    sample = 'TCTG'

    print("reference genome: ",ref_genome)
    print("sample: ",sample)
    index = kmer_Index(ref_genome,3)
    index.print()

    print("kmer matches at ")
    match = query_index(sample, ref_genome,index)
    print(match)
    
    print("exact matches at ")
    match = perfect_match(sample, ref_genome,index)
    print(match)

    
