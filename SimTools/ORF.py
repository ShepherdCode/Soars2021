class ORF ():
    def __init__(self):
        self.begin=0 # first nucleotide position, zero based
        self.end=0 # 1 + last nucleotide position
    def get_length_in_nucleotide(self):
        '''Includes the start and stop codons.'''
        return self.end-self.begin
    def get_length_in_codon(self):
        '''Includes the start and stop codons.'''
        return self.get_length_in_nucleotide()//3
    def get_begin(self):
        return self.begin
    def get_end(self):
        return self.end
