class RNA_describer():
    START='ATG'
    STOPS=['TAA','TAG','TGA']
    CLEN=3 # codon length
    def get_orf_length(self,seq):
        '''Returns non-zero only if
        given sequence begins START
        and contains a STOP at offset 3n.'''
        clen=RNA_describer.CLEN
        if len(seq)<2*clen:
            return 0 # too short
        if not seq[0:clen] == RNA_describer.START:
            return 0 # missing start
        for pos in range(0,len(seq)-clen+1,clen):
            if seq[pos:pos+clen] in RNA_describer.STOPS:
                return pos+clen
        return 0 # missing stop

    def get_longest_orf_len(self,seq):
        rlen = len(seq)
        lens=[0]*RNA_describer.CLEN
        for frame in range(0,RNA_describer.CLEN):
            start=0
        return 0

rn = RNA_describer()
print("orf len 0, no start:",rn.get_orf_length('TGAAATGA'))
print("orf len 6, start+stop:",rn.get_orf_length('ATGTGA'))
print("orf len 9, one codon:",rn.get_orf_length('ATGAAATGA'))
print("orf len 9, 3'UTR:",rn.get_orf_length('ATGAAATGAAAA'))
print("orf len 0:, no stop",rn.get_orf_length('ATGAATGA'))
#print(rn.get_longest_orf_len('ATGAAATGA'))
