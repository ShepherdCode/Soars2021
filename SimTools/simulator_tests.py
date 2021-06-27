from RNA_describe import *
from ORF_eliminator import *
from RNA_gen import *

class Simulator:
    def get_random_sequence(self):
        '''
        This is Professor Miller's class that generates a RNA sequence
        without any restriction.
        '''
        generator = Collection_Generator()
        sequence_list = generator.get_sequences() #Generates sequence in list forrmat
        sequences =''.join( sequence_list)
        #print('The sequence below is generated from Collection generator. It may or maynot have a ORF')
        return sequences
    
    def get_sequence_no_orf(self,rna):
        '''
        This function uses ORF_eliminator.
        ORF_eliminator is a class that kills any orfs created in random sequence.
        For example, 'ATGAAATAG' has a orf of length 9. My program kills the orf
        by replacing 'T' in stop codon ('TAG'in the example above) repeatedly
        unless we get orf length =0.
        There are two methods that i implemented. However, the first is more
        optimized.
        '''
        eliminator = ORF_eliminator()
        seq_with_no_orf = eliminator.eliminate_ORF(rna)
        #print(' \n The  sequence below doesnot have any orf for sure')
        return seq_with_no_orf
    
    def get_lengths_codons(self,rna):
        '''
        This functions are meant for statistical purposes.
        It will employ function which are found in RNA_describe.
        This calculates the number of start and stop codons.
        
        '''
        counter = ORF_RE()
        number_of_codons = counter.get_codon_number(rna)
        return number_of_codons
    def get_number_bases(self, rna):
        '''
         This functions are meant for statistical purposes.
        It will employ function which are found in RNA_describe.
        This calculates the number of 'A', 'C', 'T', 'G'.
        '''
        counter = ORF_RE()
        number_of_bases = counter.get_number_bases(rna)
        return number_of_bases
       
  
        
      
    
        
if __name__ == "__main__":
    simulator = Simulator()
    random_RNA = simulator.get_random_sequence()
    print('Random sequence generated from Collection_Generator')
    print(random_RNA)#Prints the random sequence generated by Collection_Genenator
    start = simulator.get_lengths_codons(random_RNA)[0]
    stop = simulator.get_lengths_codons(random_RNA)[1]
    a = simulator.get_number_bases(random_RNA)[0]
    c = simulator.get_number_bases(random_RNA)[1]
    t = simulator.get_number_bases(random_RNA)[2]
    g = simulator.get_number_bases(random_RNA)[3]


    print('THere are {} Start codons and {} Stop codons. '.format(start,stop))
    print('There are {} a, {} c, {} t, {} g ' . format(a,c,t,g))
    
    print('\n \n Sequence with no ORF from ORF_eliminator')
    
    rna = simulator.get_sequence_no_orf(random_RNA)
    start = simulator.get_lengths_codons(rna)[0]
    stop = simulator.get_lengths_codons(rna)[1]
    a = simulator.get_number_bases(rna)[0]
    c = simulator.get_number_bases(rna)[1]
    t = simulator.get_number_bases(rna)[2]
    g = simulator.get_number_bases(rna)[3]


    print('THere are {} Start codons and {} Stop codons. '.format(start,stop))
    print('There are {} a, {} c, {} t, {} g ' . format(a,c,t,g))
    print(rna)
