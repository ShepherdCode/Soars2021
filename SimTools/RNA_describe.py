class RNA_describer():
    START='ATG'
    STOPS=['TAA','TAG','TGA']
    CODON_LEN=3 # codon length
    def get_orf_length(self,seq):
        '''Count bases to end of first in-frame stop codon.
        Return count start and stop codons (min len 6).
        Assume first three bases are ATG.
        Assume we'll find a STOP codon at some 3n position.
        Else return zero.'''
        clen=RNA_describer.CODON_LEN
        if len(seq)<clen+clen:
            return 0 # too short for start + stop
        if not seq[0:clen] == RNA_describer.START:
            return 0 # missing start
        for pos in range(0,len(seq)-clen+1,clen):
            if seq[pos:pos+clen] in RNA_describer.STOPS:
                return pos+clen
        return 0 # missing stop
    def get_longest_orf(self,seq):
        '''Find longest ATG...TAG in any frame.
        Return tuple (offset,len).
        Offset starts at zero.
        Length includes ATG and TAG (min 6).
        Returns (0,0) if none found.'''
        clen=RNA_describer.CODON_LEN
        rlen = len(seq)
        longest_orf=(0,0) # offset, len
        for one_pos in range(0,rlen):
            one_len = self.get_orf_length(seq[one_pos:])
            if one_len > longest_orf[1]:
                longest_orf=(one_pos,one_len)
        return longest_orf
    def get_orf_lengths(self,seqs):
        '''Return list of largest ORF length per sequence.'''
        lens = []
        for one_seq in seqs:
            one_orf = self.get_longest_orf(one_seq)
            one_len = one_orf[1]  # orf=(offset,length)
            lens.append(one_len)
        return lens
