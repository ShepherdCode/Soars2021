import re
import math
import random
import argparse

class Random_Sequence:
    '''
    Creates the random sequence of given or default length. May/Maynot have
    ORF.
    '''
    def __init__(self):
        #Bases
        self.BASES = ['A', 'C', 'T','G']

        
    
    def generate_sequence(self,length_of_sequence,no_of_sequence):
        list = []
        counter =0
        while (counter<no_of_sequence):
            sequence = ''.join(random.choice(self.BASES) for i in range(length_of_sequence))
            list.append( sequence)
            counter+=1
        return list
    def generate_sequence_withORF(self, length_of_sequence, no_of_sequence, min_length):
        list = []
        counter =0
        orf_search = ORF_RE()
        while (counter<no_of_sequence):
            sequence = ''.join(random.choice(self.BASES) for i in range(length_of_sequence))
            if(orf_search.get_three_lengths(sequence)[1]<min_length):
                continue
            list.append( sequence)
            counter+=1
        return list
    

class ORF_RE():
    '''Regular Expressions'''
    def __init__(self):
        canonical='ATG((?!TAA|TAG|TGA)[ACGT]{3})*(TAA|TAG|TGA)'
        self.canonical_re=re.compile(canonical)
    def get_three_lengths(self,RNA):
        '''Return length of 5'UTR, longest ORF, and 3'UTR.'''
        utr5_len=0
        orf_len=0
        s=self.canonical_re.search(RNA)
        while s is not None:
            this_len=s.end()-s.start()-3 # exclude the stop codon
            if this_len>orf_len:
                orf_len=this_len
                utr5_len=s.start()
            pos = s.start()+1
            s=self.canonical_re.search(RNA,pos)
        if orf_len>0:
            utr3_len=len(RNA)-utr5_len-(orf_len+3)
        else: # no orf? then each UTR is half
            utr5_len=len(RNA)//2
            utr3_len=len(RNA)-utr5_len
        return (utr5_len,orf_len,utr3_len)
    
class ORF_eliminator:
    
        
    '''
    The following function takes RNA as a parameter
    and returns the set of lengths of 5UTF, ORF, 3UTF

    '''
    def get_coordinates(self,RNA):
        '''
        This method  is used from Professor Miller's RNA_describe.
        Professor Miller used Regular Expression to locate ORFs
        which is proven to be more effective and faster especially
        when the ORFs are in different frames.
        
        '''
        orfs = ORF_RE()
        #lengthset gives set of 3 lengths in list format.
        lengthSet = orfs.get_three_lengths(RNA)
        return (lengthSet)
   
    
    '''
------------------------------------------------------------------------------
        METHOD 1: 
        This method eliminates an ORF by killing the stop codons.
        It goes to the stop codons and replaces 'T'
        with ('A','C', 'G') randomly.
        Howver, we donot want to affect the randomness of 'T' so, whenever
        'T' is targetted, another 'T' will replace  letter in random
        place. This will make the number of 'T' constant.
------------------------------------------------------------------------------
    '''
    def eliminate_ORF(self, RNA):
        
        #choices defines the possible letters that could replace 'T'
        choices = ['A','C','G']
        lengths = self.get_coordinates(RNA) #Gets the length of 5'UTF, ORF, 3'UTF
        #lengths is the list of three length. Hence, indices are used for referral.
        utr5 = lengths[0]
        orf = lengths[1]
        utr3 = lengths[2]
        counter =0
        
        while(orf!=0):
            
            temp = list(RNA)
            temp[utr5+orf] = random.choice(choices)
            #randonNumber get the random location to inject {T}
            #This way the occurance of T isnt changed .
            
            randomNumber = random.randint(0,len(RNA)-1)
            temp[randomNumber] = 'T'
            new_RNA = "".join(temp)
            RNA = new_RNA
            lengths = self.get_coordinates(RNA)
            utr5 = lengths[0]
            orf = lengths[1]
            utr3 = lengths[2]
       # print("It took " + str(counter)+ " loops to eliminate ORF")    
        return(RNA)
        
      
        
    '''
        METHOD2:
        Injecting STOP CODONS in Middle repeatedly. 
        
    '''
        
    def eliminate_ORF2(self):
        
        
       
        orf_length = 3#randomly assigning number to get into the loop
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
        #Determines the exact position of where the STOP codon is to be injected.
            pointer = orf_starter + (pos-1)*3
            temp_RNA = temp_RNA[:pointer] + temp_codons + temp_RNA[pointer+3:]
            self.set_RNA(temp_RNA)
        
        
        return self.RNA
        
        
         
def args_parse():
        
    global args
    parser = argparse.ArgumentParser(description ='This program eliminates the orf')
    parser.add_argument('fileName',help='input filename (fasta)',type=str)
    #parser.add_argument('--resultFile', help = 'output filename (fasta)', type = str)
    args = parser.parse_args()
        
        
if __name__ == '__main__':
    
        args_parse()
        #fileName is the argument passed as a name for text file containing RNA.
        fileName = args.fileName
    
        with open(fileName, "r+") as file:
            #result.txt records the inputfile but with eliminated ORF.
            with open("result.txt", "w") as resultFile:
                lines = file.readlines()
                for line in lines:
                    if(line[0] =='>'):
                        resultFile.write(line)
                        #Breaks the Loop and continues with next line.
                        continue
                    RNA = line
                    tools = ORF_eliminator()
                    newRNA = tools.eliminate_ORF(RNA)
                    resultFile.write(newRNA)
            
        
