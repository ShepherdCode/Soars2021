import random
import string
import RNA_describe
import re

def empirical_prob_canonical_exact(orf_len):
    true=0
    count=1000000
    prob=0.0
    canonical='^ATG((?!TAA|TAG|TGA)[ACGT]{3})*(TAA|TAG|TGA)$'
    pattern=re.compile(canonical)
    letters='ACGT'
    for n in range(0,count):
        RNA=''
        for i in range(0,orf_len):
            one=random.choice(letters)
            RNA = RNA + one
        if pattern.match(RNA):
            true += 1
    prob=true/count
    print("prob = %f = %d/%d"%(prob,true,count))
    return prob

print (empirical_prob_canonical_exact(12))
