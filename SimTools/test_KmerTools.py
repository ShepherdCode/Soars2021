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
        msg="Count AA"
        assert counts['AA']==reps-1,msg
        assert freqs['AA']==pytest.approx(1.0)
        msg="Count CC"
        assert counts['CC']==0,msg
        assert freqs['CC']==pytest.approx(0.0)
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

class Test_Harvester_Counting():
    def test_harvest_1(self):
        H=True # use Harvester algorithm
        reps=10
        rna = 'A'*reps
        K = 3
        tool = KmerTools()
        empty = tool.make_dict_upto_K(K)
        counts = tool.update_count_one_K(empty,K,rna,H)
        counts = tool.harvest_counts_from_K(counts,K)
        freqs = tool.count_to_frequency(counts,K)
        msg="Count A"
        assert counts['A']==reps,msg
        assert freqs['A']==pytest.approx(1.0)
        msg="Count C"
        assert counts['C']==0,msg
        assert freqs['C']==pytest.approx(0.0)
        msg="Count AA"
        assert counts['AA']==reps-1,msg
        assert freqs['AA']==pytest.approx(1.0)
        msg="Count CC"
        assert counts['CC']==0,msg
        assert freqs['CC']==pytest.approx(0.0)
        msg="Count AAA"
        assert counts['AAA']==reps-2,msg
        assert freqs['AAA']==pytest.approx(1.0)
    def test_harvest_2(self):
        H=True # use Harvester algorithm
        rna = 'CAAAT'
        K = 3
        tool = KmerTools()
        empty = tool.make_dict_upto_K(K)
        counts = tool.update_count_one_K(empty,K,rna,H)
        counts = tool.harvest_counts_from_K(counts,K)
        freqs = tool.count_to_frequency(counts,K)
        msg="Count A"
        assert counts['A']==3,msg
        msg="Count C"
        assert counts['C']==1,msg
        msg="Count G"
        assert counts['G']==0,msg
        msg="Count T"
        assert counts['T']==1,msg
        msg="Count CA"
        assert counts['CA']==1,msg
        msg="Count AA"
        assert counts['AA']==2,msg
        msg="Count AT"
        assert counts['AT']==1,msg
        msg="Count CAA"
        assert counts['CAA']==1,msg
        msg="Count AAA"
        assert counts['AAA']==1,msg
        msg="Count AAT"
        assert counts['AAT']==1,msg
        msg="Count TAG"
        assert counts['TAG']==0,msg
