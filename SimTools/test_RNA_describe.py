import pytest
from RNA_describe import RNA_describer
from RNA_describe import ORF_counter
from RNA_describe import ORF_RE

# The following unix command will run all tests.
# $ pytest
# The -v option will list each test and show progress.
# $ pytest -v
# By default, pytest captures stdout unless the tests fail.
# Use this option to see the output of print() statements.
# $ pytest --capture=tee-sys

class Test_ORF_RE():
    def test_get_all_orfs(self):
        ore = ORF_RE()
        rna = 'CCCATGAAATGACCTGATGCCCTGACCC'
        orfs = ore.get_all_orfs(rna)
        ans = ['ATGAAATGA', 'ATGACCTGA', 'ATGCCCTGA']
        msg="Overlapping ORFs"
        assert orfs==ans,msg
    def test_get_all_orfs(self):
        ore = ORF_RE()
        rna = 'ATGCCCTGA'+'ATGCCCCCCTAG'+'CC'
        orfs = ore.get_three_lengths(rna)
        ans = (9,9,2)
        msg="Overlapping ORFs"
        assert orfs==ans,msg

class Test_ORF_counter():
    def test_three_codon_orf(self):
        oc = ORF_counter()
        msg= "Detects START CODON STOP"
        oc.set_sequence('C'+'ATG'+'CAC'+'TAG'+'C')
        assert oc.get_max_orf_len()==6,msg
        assert oc.count_maximal_orfs()==1,msg
    def test_no_codon_orf(self):
        oc = ORF_counter()
        msg = "Counts bases ATG thru TAA"
        oc.set_sequence('ATG'+'TAA'+'G')
        assert oc.get_max_orf_len()==3,msg
        oc.set_sequence('A'+'ATG'+'TAA'+'G')
        assert oc.get_max_orf_len()==3,msg
        oc.set_sequence('CA'+'ATG'+'TAA'+'G')
        assert oc.get_max_orf_len()==3,msg
    def test_no_start_codon(self):
        oc = ORF_counter()
        msg = "Detects ATG not found"
        oc.set_sequence('TGTAAGC')
        assert oc.get_max_orf_len()==0,msg
        assert oc.count_maximal_orfs()==0,msg
    def test_no_stop_codon(self):
        oc = ORF_counter()
        msg = "Detects if TAG not found"
        oc.set_sequence('ATGTACCTA')
        assert oc.get_max_orf_len()==0,msg
        assert oc.count_maximal_orfs()==0,msg
    def test_three_frames(self):
        oc = ORF_counter()
        msg = "Counts bases ATG thru TAA in frame 1"
        oc.set_sequence('CCC'+'ATG'+'AAA'+'TAA')
        assert oc.get_max_orf_len()==6,msg
        msg = "Counts bases ATG thru TAG in frame 2"
        oc.set_sequence('CC'+'ATG'+'AAA'+'TAG')
        assert oc.get_max_orf_len()==6,msg
        msg = "Counts bases ATG thru TGA in frame 3"
        oc.set_sequence('C'+'ATG'+'AAA'+'TGA')
        assert oc.get_max_orf_len()==6,msg
    def test_multiple_ORFs(self):
        oc = ORF_counter()
        msg = "Gets longest of overlapping ORFs in different frames"
        oc.set_sequence('ATG'+'AAA'+'TGA'+'AACCC'+'TGA')
        assert oc.get_max_orf_len()==9,msg
        assert oc.count_maximal_orfs()==2,msg
        msg = "Gets longest of consecutive ORFs in same frame"
        oc.set_sequence('ATG'+'TGA'+'ATG'+'AAA'+'TGA')
        assert oc.get_max_orf_len()==6,msg
        assert oc.count_maximal_orfs()==2,msg
    def test_contained_ORFs(self):
        oc = ORF_counter()
        msg = "Recognizes contained ORFs in same frame"
        oc.set_sequence('ATG'+'AAA'+'ATG'+'CCC'+'TGA')
        assert oc.get_max_orf_len()==12,msg
        assert oc.count_maximal_orfs()==1,msg
        assert oc.count_contained_orfs()==1,msg

class Test_RNA_describer():
    def test_orf_length(self):
        rn = RNA_describer()
        msg= "Require sequence starts with ATG"
        assert rn.get_orf_length('TGATGTGA')==0,msg
        msg = "Minimum requirement is start and stop"
        assert rn.get_orf_length('ATG'+'TGA')==3,msg
        msg = "Start + codon + stop = 3+3=6"
        assert rn.get_orf_length('ATG'+'AAA'+'TGA')==6,msg
        msg = "polyA tail or any 3'UTR does not count"
        assert rn.get_orf_length('ATG'+'AAA'+'TGA'+'AAAA')==6,msg
        msg = "No in-frame stop? Then no ORF"
        assert rn.get_orf_length('ATG'+'AA'+'TGA')==0,msg
    def test_longest_orf(self):
        rn = RNA_describer()
        msg = "Counts bases ATG thru TAA in frame 0"
        assert rn.get_longest_orf('ATG'+'TAA'+'G')==(0,3),msg
        msg = "Returns (0,0) if ATG not found"
        assert rn.get_longest_orf('TGTAAGC')==(0,0),msg
        msg = "Returns (0,0) if TAG not found"
        assert rn.get_longest_orf('ATGTACCTA')==(0,0),msg
        msg = "Counts bases ATG thru TAA in frame 1"
        assert rn.get_longest_orf('CCC'+'ATG'+'AAA'+'TAA')==(3,6),msg
        msg = "Counts bases ATG thru TAG in frame 2"
        assert rn.get_longest_orf('CC'+'ATG'+'AAA'+'TAG')==(2,6),msg
        msg = "Counts bases ATG thru TGA in frame 3"
        assert rn.get_longest_orf('C'+'ATG'+'AAA'+'TGA')==(1,6),msg
        msg = "Gets longest of two ORFs in same frame"
        assert rn.get_longest_orf('ATG'+'TGA'+'ATG'+'AAA'+'TGA')==(6,6),msg
        msg = "Gets longest of two ORFs in different frames"
        assert rn.get_longest_orf('ATG'+'AAA'+'TGA'+'AACCC'+'TGA')==(5,9),msg
    def test_orf_lengths(self):
        rn = RNA_describer()
        msg = "Return list of lengths"
        assert rn.get_orf_lengths(['ATG'+'TGA','ATG'+'AAA'+'TGA'])==[3,6],msg
    def test_three_lengths(self):
        rn = RNA_describer()
        msg = "ORF? Return lengths [ (5'UTR,ORF,3'UTR) ]"
        assert rn.get_three_lengths(['CAT'+'ATG'+'GGG'+'TGA'+'AAA'])==[(3,6,3)],msg
        msg = "No ORF? Return lengths [ (half,0,half) ]"
        assert rn.get_three_lengths(['CCC'+'AAA'])==[(3,0,3)],msg
