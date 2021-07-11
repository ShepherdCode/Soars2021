#!/usr/bin/env python
# coding: utf-8

# # MLP GenCode 
# Wen et al 2019 used DNN to distinguish GenCode mRNA/lncRNA.
# Based on K-mer frequencies, K={1,2,3}, they reported 99% accuracy.
# Their CNN used 2 Conv2D layers of 32 filters of width 3x3, max pool 2x2, 25% drop, dense 128.
# Can we reproduce that with MLP layers instead of CNN?
# Extract features as list of K-mer frequencies for K={1,2,3}.

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
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/GenCodeTools.py')
    with open('GenCodeTools.py', 'w') as f:
        f.write(r.text)  
    from GenCodeTools import GenCodeLoader
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/KmerTools.py')
    with open('KmerTools.py', 'w') as f:
        f.write(r.text)  
    from KmerTools import KmerTools
else:
        print("CoLab not working. On my PC, use relative paths.")
        DATAPATH='data/'  # must end in "/"
        sys.path.append("..") # append parent dir in order to use sibling dirs
        from SimTools.GenCodeTools import GenCodeLoader
        from SimTools.KmerTools import KmerTools
MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login


# ## Data Load
# Restrict mRNA to those transcripts with a recognized ORF.

# In[5]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'
PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME


# In[6]:


# Full GenCode ver 38 human is 106143 pc + 48752 nc and loads in 7 sec.
# Expect fewer transcripts if special filtering is used.
loader=GenCodeLoader()
loader.set_label(1)
loader.set_check_utr(True)
pcdf=loader.load_file(PC_FULLPATH)
print("PC seqs loaded:",len(pcdf))
loader.set_label(0)
loader.set_check_utr(False)
ncdf=loader.load_file(NC_FULLPATH)
print("NC seqs loaded:",len(ncdf))
show_time()


# ## Data Prep

# In[7]:


# Length filter
def apply_length_filter(df,low,high):
    # The pandas query language is strange, 
    # but this is MUCH faster than loop & drop.
    return df[ (df['seqlen']>=low) & (df['seqlen']<=high) ]

PC_LENGTHS=(200,4000)
NC_LENGTHS=(200,4000)  # oddly, the paper used 250 to 3500 for NC only
pc_filtered=apply_length_filter(pcdf,PC_LENGTHS[0],PC_LENGTHS[1])
nc_filtered=apply_length_filter(ncdf,NC_LENGTHS[0],NC_LENGTHS[1])
show_time()
print("PC seqs pass filter:",len(pc_filtered),"max len:",pc_filtered['seqlen'].max())
print("NC seqs pass filter:",len(nc_filtered),"max len:",nc_filtered['seqlen'].max())
# Garbage collection to reduce RAM footprint
pcdf=None
ncdf=None


# In[9]:


# Random selections of same number of sequences
def shuffle_rows(df,and_index=False):
    if and_index:
        # This is the default and the only choice in old Pandas.
        # The index gets shuffled too, so df.iloc[0] has index != 0. Strange!
        return df.sample(frac=1)
    # This is new in Pandas 1.3.
    # After shuffling, df.iloc[0] has index == 0.
    return df.sample(frac=1,ignore_index=True)
pc_filtered = suffle_rows(pc_filtered)
pc_filtered = suffle_rows(pc_filtered)


# In[10]:


MAX_K = 3
tool = KmerTools()
pc_counts = tool.make_dict_upto_K(MAX_K)
for sample in pc_all:
    print(sample)
    tool.update_count_one_K(pc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(pc_counts,MAX_K)
print("PC counts:\n",pc_counts)
pc_freqs = tool.count_to_frequency(pc_counts,MAX_K)
print ("Frequency:\n",pc_freqs)


# In[ ]:


nc_counts = tool.make_dict_upto_K(MAX_K)
for sample in nc_all:
    tool.update_count_one_K(nc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(nc_counts,MAX_K)
print("NC counts:\n",nc_counts)
nc_freqs = tool.count_to_frequency(nc_counts,MAX_K)
print ("Frequency:\n",nc_freqs)


# In[ ]:


X,y = prepare_inputs_len_x_alphabet(pc_train,nc_train,ALPHABET) # shuffles
print("Data ready.")

