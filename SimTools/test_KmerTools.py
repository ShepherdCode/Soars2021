import pytest
from KmerTools import KmerTools

# The following unix command will run all tests.
# $ pytest
# The -v option will list each test and show progress.
# $ pytest -v
# By default, pytest captures stdout unless the tests fail.
# Use this option to see the output of print() statements.
# $ pytest --capture=tee-sys

class Test_Naive_Counting():
    def test_naive_1(self):
        H=False # avoid Harvester algorithm
        reps=10
        rna = 'A'*reps
        K = 3
        tool = KmerTools()
        empty = tool.make_dict_upto_K(K)
        counts = tool.update_count_upto_K(empty,K,rna,H)
        freqs = tool.count_to_frequency(counts,K)
        msg="Count A"
        assert counts['A']==reps,msg
        assert freqs['A']==pytest.approx(1.0)
        msg="Count C"
        assert counts['C']==0,msg
        assert freqs['C']==pytest.approx(0.0)
    def test_naive_2(self):
        H=False # avoid Harvester algorithm
        reps=5
        rna = 'AT'*reps
        K = 3
        tool = KmerTools()
        empty = tool.make_dict_upto_K(K)
        counts = tool.update_count_upto_K(empty,K,rna,H)
        freqs = tool.count_to_frequency(counts,K)
        msg="Count AT"
        assert counts['AT']==reps,msg
        assert freqs['AT']==pytest.approx(5.0/9.0)
        msg="Count TA"
        assert counts['TA']==reps-1,msg
        assert freqs['TA']==pytest.approx(4.0/9.0)
        msg="Count AA"
        assert counts['AA']==0,msg
        assert freqs['AA']==pytest.approx(0.0)
        msg="Count CG"
        assert counts['CG']==0,msg
        assert freqs['CG']==pytest.approx(0.0)
    def test_naive_3(self):
        H=False # avoid Harvester algorithm
        reps=4
        rna = 'ATG'*reps
        K = 3
        tool = KmerTools()
        empty = tool.make_dict_upto_K(K)
        counts = tool.update_count_upto_K(empty,K,rna,H)
        freqs = tool.count_to_frequency(counts,K)
        msg="Count ATG"
        assert counts['ATG']==reps,msg
        msg="Count AT"
        assert counts['AT']==reps,msg
        msg="Count TG"
        assert counts['GA']==reps-1,msg
        msg="Count TA"
        assert counts['TA']==0,msg
