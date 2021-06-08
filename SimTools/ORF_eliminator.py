from RNA_describe import RNA_describer
from RNA_describe import ORF_counter
from RNA_describe import ORF_RE
import math
import random
'''
This class constructs the RNA sequence without 

'''
class ORF_eliminator:
    def __init__(self, RNA):
        self.RNA =RNA
    def print_RNA(self):
        print(self.RNA)
    def get_coordinates(self):
        orfs = ORF_RE()
        lengthSet = orfs.get_three_lengths(self.RNA)
        return (lengthSet)
    def set_RNA(self, newRNA):
        self.RNA = newRNA
    
    '''
        METHOD 1: 
        This method eliminates an ORF by killing the stop codons.
        It goes to the stop codons and replaces 'T' with ('A','C', 'G') randomly
    '''
    def eliminate_ORF(self):
        #For replacing 'T' in stop codons.
        
        choices = ['A','C','G']
        lengths = self.get_coordinates()
        utr5 = lengths[0]
        ofr = lengths[1]
        utr3 = lengths[2]
        #Repeats the process unless the ORF comes down to 0.
        while(ofr!=0):
            temp = list(self.RNA)
            temp[utr5+ofr] = random.choice(choices)
            new_RNA = "".join(temp)
            print(new_RNA)
            self.set_RNA(new_RNA)
            lengths = self.get_coordinates()
            utr5 = lengths[0]
            ofr = lengths[1]
            utr3 = lengths[2]
            
        return(self.RNA)
        
        
        
        '''
        METHOD2:
        Injecting STOP CODONS in Middle repeatedly. 
        
        '''
        
    def eliminate_ORF2(self):
        
        lengths= self.get_coordinates();
        #return lengths
        #Working for coordinates
        orf_length = 3
        while(orf_length!=0):
            lengths= self.get_coordinates();
            orf_length = lengths[1]
            orf_starter = lengths[0]
            orf_end = lengths[0]+lengths[1]+3
            no_of_codons = int((orf_end-orf_starter)/3)
        
            STOP_CODONS = ['TAG','TGA','TAA']
        #Determines mid-codon. 
            pos= math.ceil(no_of_codons/2)
            temp_codons = random.choice(STOP_CODONS)
            temp_RNA = self.RNA
        #Determines the exact position of codons.
            pointer = orf_starter + (pos-1)*3
            temp_RNA = temp_RNA[:pointer] + temp_codons + temp_RNA[pointer+3:]
            self.set_RNA(temp_RNA)
        
        
        return self.RNA
        
        
         
        
        
if __name__ == '__main__':
     #RNA = 'TATAATGTCACTGTTTTACTATAGGGTCGTAGGGTACTGGACCCCGCATT'
     #RNA = 'ATGTTGAGCCCCATTAGATTACCAGATTCTGACCGGATCCGTTTAATAGA'
     RNA = 'AGGATGCCCTGA'+'ATGCCCCCCTAG'+'CC'
     trial = ORF_eliminator(RNA)
  
     print(trial.eliminate_ORF2())
   