#!/usr/bin/env python
# coding: utf-8

# # Probability of an ORF
# 
# For getting started, we'll solve a simple version.
# Assume the RNA length is a multiple of 3.
# Only count ORFs in frame 1.
# Assume every codon is a random uniform selection from 64 choices.  
# Solve one one exact length in codons.
# Include every ORF of length >= codons/2. 

# In[68]:


import time
def show_time():
    t = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t)))
show_time()


# Let  
# W = P(NonStart) = 63/64  
# M = P(Start) = 1/64  
# S = P(Stop) = 3/64  
# A = P(NonStop)= 61/64  
# C = P(AnyCodon) = 64/64  
# n = Length of RNA in bases.  
# L = ORF length in codons (includes start, excludes stop)  
# H = ceil(n/2)  
# I = floor(n/2)  
# Porf(n) = Probability of an ORF, L>=n/2.  
# 
# $P_{orf}(n)=\sum_{i=1}^I\sum_{j=i}^L[W^{i-1}MA^{j}S]$

# In[69]:


import math
P_NoStart = 63/64
P_Start = 1/64
P_Stop = 3/64
P_Amino = 61/64
P_Codon = 64/64


# In[70]:


def get_prob_by_sum(rna_len):
    psum = 0.0  # sum of probs
    min_cds = math.ceil(rna_len/2)
    min_amino = min_cds-2
    max_start = math.floor(rna_len/2)
    for start_pos in range(1,max_start+1):
        for stop_pos in range(start_pos+min_cds,rna_len+1):
            num_pre = start_pos-1
            num_amino = stop_pos-start_pos-1
            pone = (P_NoStart**num_pre)*P_Start*(P_Amino**num_amino)*P_Stop
            psum += pone
    return psum


# In[78]:


print("  RNA  CODONS  PROB")
for c in range(1,11,1):
    rna_codons=c
    rna_len=rna_codons*3
    porf = get_prob_by_sum(rna_codons)
    print("%5d  %6d  %.5f"%(rna_len,rna_codons,porf))


# In[79]:


print("  RNA  CODONS  PROB")
for c in range(33,700,33):
    rna_codons=c
    rna_len=rna_codons*3
    porf = get_prob_by_sum(rna_codons)
    print("%5d  %6d  %.5e"%(rna_len,rna_codons,porf))


# In[ ]:




