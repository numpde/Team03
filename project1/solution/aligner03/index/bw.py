# Can be used to generated the Burrows Wheeler transformation (BWT) and the corresponding offsets in the suffix array (OSA). Upon creation of the object the BWT and OSA are created and stored (Both of size len(reference_genome) ). A method to decode a BWT is provided. Works with strings of any shape. The encoding can be done via suffix arrays (default) or via rotation of the string (specified at initialisation).

# no _
class Burrows_Wheeler(object):

    def __init__(self, reference_genome, strategy='S'):
        # choose encoding strategy to build the Burrows-Wheeler Transform
        # R = via rotation of the string
        # S = via suffix array (preferred method)

        if strategy not in ['S', 'R']:
            raise ValueError('strategy needs to be S (suffix array) or R (rotation)')

        if len(reference_genome) < 1:
            raise ValueError('please provide a non empty string for the reference genome')

        self.encode_strategy = strategy

        # get Burrows Wheeler transformation and corresponding offsets in the suffix array
        (self.code, self.sa) = self.bwt_encode(reference_genome)

    # Given a reference genome, returns the Burrows Wheeler transformation
    def bwt_encode(self, reference_genome):

        if (self.encode_strategy == 'R'):
            return self.bwt_encode_rot(reference_genome)
        else:
            return self.bwt_encode_suffix(reference_genome)

    # Given a reference_genome, returns the Burrows-Wheeler Transform constructed by the rotation of the reference_genome
    def bwt_encode_rot(self, reference_genome):

        n = len(reference_genome)
        bw_transform = sorted([reference_genome[i:n] + reference_genome[0:i] for i in range(n)])
        L = ''.join(last[-1] for last in bw_transform)

        return L

        # Given a reference_genome, returns the Burrows-Wheeler Transform constructed by suffix array offsets

    def bwt_encode_suffix(self, reference_genome):

        suffix_array_offsets = self.suffix_array_offsets(reference_genome)
        bw_transform = []

        for w in suffix_array_offsets:
            if w == 0:
                bw_transform.append('$')
            else:
                bw_transform.append(reference_genome[w - 1])

        return (''.join(bw_transform), suffix_array_offsets)

    # Given a reference_genome, returns the sorted suffix array (suffix, offset)    
    def suffix_array(self, reference_genome):

        n = len(reference_genome)
        suffix_array = sorted([(reference_genome[i:], i) for i in range(n)])

        return suffix_array

    # Given a reference_genome, returns the offsets of the sorted suffixes in reference genome
    def suffix_array_offsets(self, reference_genome):

        suffix_array = self.suffix_array(reference_genome)
        offsets = [i[1] for i in suffix_array]

        return (offsets)

    # returns the original string
    def decode(self, bw_transform=None):

        bw_transform = bw_transform or self.code

        matrix = [''] * len(bw_transform)

        # add Burrows Wheeler transform as first column and sort
        for i in range(len(bw_transform)):
            matrix = sorted(bw_transform[i] + matrix[i] for i in range(len(bw_transform)))

        s = ''
        for row in matrix:
            if row.endswith("$"):
                s = row

        return s


if __name__ == "__main__":
    ref_genome = 'This is a test string to verify the correctness of the Burrows Wheeler transformation. This functionality is used for the 1. project of the course Computational Biomedicine at ETH Zurich.'

    bwt = Burrows_Wheeler(ref_genome)

    print("test: \n", ref_genome, "\n")
    print("encoded: \n", bwt.code, "\n")
    print("decoded:\n", bwt.decode(), "\n")
