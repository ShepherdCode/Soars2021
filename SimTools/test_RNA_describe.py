import pytest
from RNA_describe import RNA_describer

class Test_RNA_describer():
    def test_orf_length(self):
        rn = RNA_describer()
        # Require sequence starts with ATG
        assert(rn.get_orf_length('TGATGTGA')==0)
        # Minimum requirement is start and stop
        assert(rn.get_orf_length('ATG'+'TGA')==6)
        # Start + codon + stop = 3*3=9
        assert(rn.get_orf_length('ATG'+'AAA'+'TGA')==9)
        # polyA tail or any 3'UTR does not count
        assert(rn.get_orf_length('ATG'+'AAA'+'TGA'+'AAAA')==9)
        # No in-frame stop? Then no ORF
        assert(rn.get_orf_length('ATG'+'AA'+'TGA')==0)

#assert(rn.get_longest_orf_len('ATGAAATGA'))
