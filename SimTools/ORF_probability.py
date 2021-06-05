import random
import string
import RNA_describe
import re

class ORF_probability():
    def canonical_exact(self,rna_bases):
        prob=0.0
        if rna_bases%3==0 and rna_bases>=6:
            prob = 1/64 * 3/64 # start and stop codons
            pResidue = 61/64 # prob of non-stop codon
            if rna_bases>6:
                codons=(rna_bases-6)//3
                pns=pResidue**codons
                prob=prob*pns
        return prob
    def canonical_range(self,min_orf_len,rna_len):
        prob=0.0
        sum=0.0
        if rna_len>=min_orf_len+3:
            prob = 1/64 * 3/64 # start and stop codons
            pResidue = 61/64 # prob of non-stop codon
            n=(rna_len-3)//3
            m=min_orf_len//3
            print("m,n=",m,n)
            for i in range(m,n+1):
                placements=n-i+1
                pns=pResidue**(i-1)
                sum = sum+placements*pns
        prob = prob * sum
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
