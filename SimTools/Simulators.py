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
    
if __name__ == '__main__':
    print('Working')
    PC=True
    NC=False
    RNA_LEN=150
    ks = KeagenSimulator()
    seq1 = ks.generate_seq(RNA_LEN, PC,3)
    print(seq1)
    
    