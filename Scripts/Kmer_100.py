#!/usr/bin/env python
# coding: utf-8

# # K-mer 
# Basic K-mer counting.

# In[1]:


import time
def show_time():
    t = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t)))
show_time()


# In[2]:


PC_SEQUENCES=32000
NC_SEQUENCES=32000
RNA_LEN=32
CDS_LEN=16


# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[4]:


import sys
IN_COLAB = False
try:
    from google.colab import drive
    IN_COLAB = True
except:
    pass
if IN_COLAB:
    print("On Google CoLab, mount cloud-local file, get our code from GitHub.")
    PATH='/content/drive/'
    #drive.mount(PATH,force_remount=True)  # hardly ever need this
    #drive.mount(PATH)    # Google will require login credentials
    DATAPATH=PATH+'My Drive/data/'  # must end in "/"
    import requests
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_describe.py')
    with open('RNA_describe.py', 'w') as f:
        f.write(r.text)  
    from RNA_describe import Random_Base_Oracle
else:
        print("CoLab not working. On my PC, use relative paths.")
        DATAPATH='data/'  # must end in "/"
        sys.path.append("..") # append parent dir in order to use sibling dirs
        from SimTools.RNA_describe import Random_Base_Oracle
MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login


# ## K-mer counting

# ### Functions to create the dict of {kmer:count}

# In[5]:


def make_kmer_keys(K):
    shorter_kmers=['']
    for i in range(K):
        longer_kmers=[]
        for mer in shorter_kmers:
            # No support for N or any non-ACGT bases.
            longer_kmers.append(mer+'A')
            longer_kmers.append(mer+'C')
            longer_kmers.append(mer+'G')
            longer_kmers.append(mer+'T')
        shorter_kmers = longer_kmers
    return shorter_kmers
def make_kmer_dict(keys,init=0):
    return dict.fromkeys(keys,init)
def make_dict_upto_K(max_K):
    keys=make_kmer_keys(1)
    for k in range(2,max_K+1):
        keys.extend(make_kmer_keys(k))
    counts = make_kmer_dict(keys)
    return counts


# ### Naive K-mer counting algorithm
# Algorithm:  
# 1. for every string  
#     1. for every K  
#         1. for every position  
#             1. kmer=substring
#             2. count{kmer}++

# In[6]:


def update_count_one_K(counts,K,rna,tail=False):
    L = len(rna)
    padding=" "*(K-1)
    padded=rna+padding
    for i in range(0,L-K+1):
        kmer=padded[i:i+K]
        counts[kmer] += 1
    if tail and K>1:  
        # for Harvester algorithm, count last letters as special case
        for start_pos in range(L-K+1,L):
            for end_pos in range(start_pos+1,L+1):
                kmer=rna[start_pos:end_pos]
                counts[kmer] += 1
    return counts
def update_count_upto_K(counts,max_K,sample,tail=False):
    for i in range(1,max_K+1):
        update_count_one_K(counts,i,sample,tail)
    return counts


# ### Harvester K-mer counting algorithm
# Algorithm:  
# 1. Count K-mers for max K only  
# 2. For each K-mer in counts table:  
#     1. For every prefix of the K-mer:  
#         1. count{prefix} += count{kmer}  
# 3. Handle last K-1 letters of each string as special case

# In[7]:


def harvest_counts_from_K(counts,max_K):
    for kmer in counts.keys():
        klen = len(kmer)
        kcnt = counts[kmer]
        if klen==max_K and kcnt>0:
            for i in range(1,klen):
                prefix = kmer[:i]
                counts[prefix] += kcnt
    return counts


# In[8]:


def count_to_frequency(counts,max_K):
    freqs = dict.fromkeys(counts.keys(),0.0)
    for k in range(1,max_K+1):
        tot = 0
        for kmer in counts.keys():
            if len(kmer)==k:
                tot += counts[kmer]
        for kmer in counts.keys():
            if len(kmer)==k:
                freqs[kmer] = 1.0*counts[kmer]/tot
    return freqs


# ## Demo

# ### Demo: Naive algorithm

# In[9]:


MAX_K = 3
counts1 = make_dict_upto_K(MAX_K)
print("Initial counts:\n",counts1)

sample = "ACCGGGTTTTACGTACGT"
update_count_upto_K(counts1,MAX_K,sample)
print("Final counts:\n",counts1)


# ### Demo: Harvester algorithm

# In[10]:


MAX_K = 3
counts2 = make_dict_upto_K(MAX_K)
print("Initial counts:\n",counts2)

sample = "ACCGGGTTTTACGTACGT"
update_count_one_K(counts2,MAX_K,sample,True)
print("Partial counts (just max K and special case letters)\n:",counts2)
harvest_counts_from_K(counts2,MAX_K)
print("Final counts (includes smaller values of K):\n",counts2)


# In[11]:


if counts1==counts2:
    print("Success. Harvester output matches naive results!")
else:
    print("Fail. Harvester output differs from naive results!")


# In[12]:


freqs = count_to_frequency(counts2,MAX_K)
print ("Frequency:\n",freqs)


# ## Demo on large dataset

# In[13]:


rbo=Random_Base_Oracle(RNA_LEN,True)
pc_all,nc_all = rbo.get_partitioned_sequences(CDS_LEN,10) # just testing
pc_all,nc_all = rbo.get_partitioned_sequences(CDS_LEN,PC_SEQUENCES)
print("Use",len(pc_all),"PC seqs")
print("Use",len(nc_all),"NC seqs")


# In[14]:


MAX_K = 3
pc_counts = make_dict_upto_K(MAX_K)
for sample in pc_all:
    update_count_one_K(pc_counts,MAX_K,sample,True)
harvest_counts_from_K(pc_counts,MAX_K)
print("PC counts:\n",pc_counts)
pc_freqs = count_to_frequency(pc_counts,MAX_K)
print ("Frequency:\n",pc_freqs)


# In[15]:


nc_counts = make_dict_upto_K(MAX_K)
for sample in nc_all:
    update_count_one_K(nc_counts,MAX_K,sample,True)
harvest_counts_from_K(nc_counts,MAX_K)
print("NC counts:\n",nc_counts)
nc_freqs = count_to_frequency(nc_counts,MAX_K)
print ("Frequency:\n",nc_freqs)


# In[ ]:




