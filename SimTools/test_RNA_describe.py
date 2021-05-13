import pytest
from RNA_describe import RNA_describer

class Test_RNA_describer():
    def test_orf_length(self):
        rn = RNA_describer()
        msg= "Require sequence starts with ATG"
        assert rn.get_orf_length('TGATGTGA')==0,msg
        msg = "Minimum requirement is start and stop"
        assert rn.get_orf_length('ATG'+'TGA')==6,msg
        msg = "Start + codon + stop = 3*3=9"
        assert rn.get_orf_length('ATG'+'AAA'+'TGA')==9,msg
        msg = "polyA tail or any 3'UTR does not count"
        assert rn.get_orf_length('ATG'+'AAA'+'TGA'+'AAAA')==9,msg
        msg = "No in-frame stop? Then no ORF."
        assert rn.get_orf_length('ATG'+'AA'+'TGA')==0,msg
    def test_longest_orf(self):
        rn = RNA_describer()
        msg = "Counts bases ATG thru TAA in frame 0"
        assert rn.get_longest_orf('ATG'+'TAA'+'G')==(0,6),msg
        msg = "Returns (0,0) if ATG not found."
        assert rn.get_longest_orf('TGTAAGC')==(0,0),msg
        msg = "Returns (0,0) if TAG not found."
        assert rn.get_longest_orf('ATGTACCTA')==(0,0),msg
        msg = "Counts bases ATG thru TAA in frame 1"
        assert rn.get_longest_orf('CCC'+'ATG'+'AAA'+'TAA')==(3,9),msg
        msg = "Counts bases ATG thru TAG in frame 2"
        assert rn.get_longest_orf('CC'+'ATG'+'AAA'+'TAG')==(2,9),msg
        msg = "Counts bases ATG thru TGA in frame 3"
        assert rn.get_longest_orf('C'+'ATG'+'AAA'+'TGA')==(1,9),msg

#assert(rn.get_longest_orf_len('ATGAAATGA'))
