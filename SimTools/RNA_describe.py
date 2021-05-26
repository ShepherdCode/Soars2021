def assert_imported_RNA_describe():
    return True

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
        '''Given list of zero to many sequences.
        Return list of largest ORF length per sequence.'''
        lens = []
        for one_seq in seqs:
            # orf=(offset,length)
            one_orf = self.get_longest_orf(one_seq)
            one_len = one_orf[1]
            lens.append(one_len)
        return lens
    def get_three_lengths(self,seqs):
        '''Given list of zero to many sequences.
        Find and use the longest ORF per sequence.
        Return list of three lengths per sequence, like this:
        [ (5'UTR, ORF, 3'UTR) ].
        This support statistical characterization of a data set.
        For sequences lacking an ORF, return (len/2,0,len/2).'''
        bounds = []
        for one_seq in seqs:
            # orf=(offset,length)
            one_orf = self.get_longest_orf(one_seq)
            seq_length=len(one_seq)
            orf_length=one_orf[1]
            if orf_length==0:
                one_bound = ((seq_length+1)//2,0,seq_length//2)
            else:
                utr5_length=one_orf[0]
                utr3_length=seq_length-utr5_length-orf_length
                one_bound = (utr5_length,orf_length,utr3_length)
            bounds.append(one_bound)
        return bounds
