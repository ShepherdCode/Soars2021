import pytest
from ORF_probability import *

class Test_ORF_probability():
    def in_range(self,number,low,high):
        return number>=low and number<=high
    def test_canonical_exact(self):
        op=ORF_probability()
        msg="p(ORF in 6 bases)"
        prob=op.canonical_exact(6)
        assert self.in_range(prob,0.0007,.0008),msg
        msg="p(ORF in 12 bases) analytical"
        prob=op.canonical_exact(12)
        assert self.in_range(prob,0.0006,.0007),msg
        msg="p(ORF in 12 bases) by random choice - results may vary"
        prob=op.empirical_prob_canonical_exact(12,500000)
        assert self.in_range(prob,0.0006,.0007),msg
