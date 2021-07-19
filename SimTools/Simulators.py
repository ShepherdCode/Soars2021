import random
import itertools

class KeagenSimulator():
    def __init__(self):
        self.Pairs = ['A', 'C', 'T', 'G']
        self.Codons = itertools.product(self.Pairs, repeat=3)
        self.all_codons = ["".join(codon) for codon in self.Codons]
        self.start = "ATG"
        self.stop = ("TAA", "TAG", "TGA") 
        self.no_stops = [codon for codon in self.all_codons if codon not in self.stop]
    
    def shift_frame(self,input_seq,frame = 2):
      output = input_seq
      if frame in (1,2,3):
        for i in range(frame-1):
          output.insert(0, random.choice(("A","G","C","T")))
          return output
      else:
        raise ValueError("Frame Must Be 1, 2 or 3. Frame Entered: " +frame)
    def codon_placer(self,length, coding = True):
      lno_stops = self.no_stops
      lall_codons = self.all_codons

      if coding == True:
        return random.choices(lno_stops, k=length)
      else:
        return random.choices(lall_codons,k=length)
    
    def get_index_placement(self,total_codons):
      quintile = total_codons // 5
      variation = quintile // 3
      start_placement = quintile + random.randint(-variation, variation)
      stop_placement = (quintile*4) + random.randint(-variation, variation)
      return start_placement, stop_placement
    
    def generate_seq(self,length, coding = False, frame = 1):
        codons_to_place = (length//3) + 1

        if coding and frame in (1,2,3):
            start_index, stop_index = self.get_index_placement(codons_to_place)
            UTR_5_len = start_index-1
            orf_length = stop_index-start_index - 2
            UTR_3_len = codons_to_place - stop_index + 1
            UTR_5 = self.codon_placer(UTR_5_len, False)
            sequence_orf = self.codon_placer(orf_length, True)
            sequence_orf.insert(0, self.start)
            sequence_orf.append(random.choice(self.stop))
            UTR_3 = self.codon_placer(UTR_3_len, False)

            UTR_5.extend(sequence_orf)
            UTR_5.extend(UTR_3)
            output = self.shift_frame(UTR_5, frame)
            output = ''.join(output)
            # Unclear what last return pair represented.
            #return output[0:length], coding, frame, (start_index, stop_index)
            return output[0:length], coding, frame

        elif not coding and frame in (1,2,3):
          placed_codons = codon_placer(codons_to_place, coding)
          output = shift_frame(placed_codons, frame)
          output = ''.join(output)
          return output[0:length] , coding, frame
        else:
          raise ValueError("Frame must be 1, 2 or 3")
        
class SamonSimulator():
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
    def eliminate_ORF(self, RNA,CDS_length):
        
        #choices defines the possible letters that could replace 'T'
        choices = ['A','C','G']
        lengths = self.get_coordinates(RNA) #Gets the length of 5'UTF, ORF, 3'UTF
        #lengths is the list of three length. Hence, indices are used for referral.
        utr5 = lengths[0]
        orf = lengths[1]
        utr3 = lengths[2]
        counter =0
        
        while(orf>CDS_length):
            
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
        
    def eliminate_ORF2(self,RNA,CDS_LENGTH):
        
        
       
        orf_length = 3#randomly assigning number to get into the loop
        while(orf_length>CDS_LENGTH):
            lengths= self.get_coordinates(RNA);
            orf_length = lengths[1]
            orf_starter = lengths[0]
            orf_end = lengths[0]+lengths[1]+3
            no_of_codons = int((orf_end-orf_starter)/3)
        
            STOP_CODONS = ['TAG','TGA','TAA']
        #Determines mid-codon. 
            pos= math.ceil(no_of_codons/2)
            temp_codons = random.choice(STOP_CODONS)
            temp_RNA = RNA
        #Determines the exact position of where the STOP codon is to be injected.
            pointer = orf_starter + (pos-1)*3
            temp_RNA = temp_RNA[:pointer] + temp_codons + temp_RNA[pointer+3:]
            RNA=temp_RNA
        
        
        return RNA
    
if __name__ == '__main__':
    print('Working')
    PC=True
    NC=False
    RNA_LEN=150
    ks = KeagenSimulator()
    seq1 = ks.generate_seq(RNA_LEN, PC,3)
    print(seq1)
    
    