import random
import string
import RNA_describe
import re

class ORF_probability():
    def canonical_exact(self,orf_len):
        prob=0.0
        pStop = 3/64 # prob of stop codon
        pStart = 1/64 # prob of start codon
        pNonStop = 61/64 # prob of non-stop codon
        if orf_len%3==0:
            if orf_len>=3:
                prob=pStart*pStop
                if orf_len>3:
                    codons=(orf_len-6)//3
                    pns=pNonStop**codons
                    prob=prob*pns
        return prob
    def canonical(self):
        for stop_pos in range(min_orf_len,seq_len-3+1):
            print("stop pos=",stop_pos)
            for start_pos in range(stop_pos-3,-1,-3):
                print("start pos=",stop_pos)
                orf_len=stop_pos-start_pos
                if (orf_len >= min_orf_len):
                    p1=pStart*pStop
                    prob += p1
        print(prob)
        return prob
    def start_codon(self,given_length):
        return 0.0
    def no_stop_codon(self,given_length):
        return (3/64)**given_length
    def empirical_prob_canonical_exact(self,orf_len,trials=1000):
        '''By random trials, empirically compute
        the probability of a canonical ORF
        that is exactly the given length.'''
        true=0
        prob=0.0
        canonical='^ATG((?!TAA|TAG|TGA)[ACGT]{3})*(TAA|TAG|TGA)$'
        pattern=re.compile(canonical)
        letters='ACGT'
        for n in range(0,trials):
            RNA=''
            for i in range(0,orf_len):
                one=random.choice(letters)
                RNA = RNA + one
            if pattern.match(RNA):
                true += 1
        prob=true/trials
        print("prob = %f = %d/%d"%(prob,true,trials))
        return prob
