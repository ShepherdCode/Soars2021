def assert_imported_RNA_describe():
    return True

class ORF_counter():
    '''Assume RNA composed of ACGT upper case no N.'''
    def __init__(self,seq):
        self.RNA=seq
        self.max_orf_len=0
        self.num_maximal_orfs=0
        self.num_contained_orfs=0
        self.prev_end=[0,0,0]
        self.prev_start=[0,0,0]
        self.__update__()
    def get_max_orf_len(self):
        return self.max_orf_len
    def count_maximal_orfs(self):
        return self.num_maximal_orfs
    def count_contained_orfs(self):
        return self.num_contained_orfs
    def __update__(self):
        RNA=self.RNA
        print("RNA is",RNA)
        pos=len(RNA)-3
        if pos<3:
            return  # smallest ORF is 6 e.g. ATG TAA
        AG = ['A','G']
        while(pos>0): # stop one before first letter
            minus1=RNA[pos-1]
            base=RNA[pos]
            plus1=RNA[pos+1]
            plus2=RNA[pos+2]
            print("base at",pos,"are",minus1,base,plus1,plus2)
            if base=='T' and plus1 in AG and plus2 in AG:
                if plus1=='G':
                    if plus2=='A':
                        self.__orf_ends_at__(pos) # TGA
                    if minus1=='A':
                        self.__orf_starts_at__(pos-1) # ATGX
                else: #   TAA, TAG
                    self.__orf_ends_at__(pos)
                pos=pos-3
            elif base=='C':
                pos=pos-3
            else:
                if base=='A' and plus1=='T' and plus2=='G':
                    self.__orf_starts_at__(pos) # ATGX
                pos=pos-1
    def __orf_ends_at__(self,pos):
        # to do: keep track of frame to reduce modulus calls
        frame=pos%3
        print("stop at",pos,"frame",frame)
        self.prev_start[frame]=0
        self.prev_end[frame]=pos
    def __orf_starts_at__(self,pos):
        # to do: keep track of frame to reduce modulus calls
        frame=pos%3
        print("start at",pos,"frame",frame)
        prev_start=self.prev_start[frame]
        prev_end=self.prev_end[frame]
        if prev_end>0:
            if prev_start>0:
                # previously reported ORF wasn't maximal afterall
                self.num_contained_orfs += 1
            else:
                # this orf is maximal (for now)
                self.num_maximal_orfs += 1
            this_len=prev_end-pos+3
            if this_len>self.max_orf_len:
                self.max_orf_len=this_len
        self.prev_start[frame]=pos

class RNA_describer():
    '''Assume RNA composed of ACGT upper case no N.'''
    START='ATG'
    STOPS=['TAA','TAG','TGA']
    CODON_LEN=3 # codon length
    def get_orf_length(self,seq):
        '''Count bases to end of first in-frame stop codon.
        Return count start and stop codons (min len 6).
        Return 0 unless first three bases are ATG.
        Assume we'll find a STOP codon at some 3n position.
        Return 0 if no in-frame stop is found.'''
        clen=RNA_describer.CODON_LEN
        if seq[0:clen] == RNA_describer.START:
            if len(seq)>=clen+clen:
                for pos in range(clen,len(seq)-clen+1,clen):
                    if seq[pos:pos+clen] in RNA_describer.STOPS:
                        return pos+clen
        return 0
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
            # TO DO: an optimization for seq like START,START,STOP.
            # No need to count length at second START.
            if seq[one_pos:one_pos+clen] == RNA_describer.START:
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
