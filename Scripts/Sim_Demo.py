#!/usr/bin/env python
# coding: utf-8

# # Sim_Demo
# Demonstrate the RNA simulator.

# In[1]:


import sys
try:
    from google.colab import drive
    IN_COLAB = True
    print("On Google CoLab, mount cloud-local file, get our code from GitHub.")
    PATH='/content/drive/'
    #drive.mount(PATH,force_remount=True)  # hardly ever need this
    #drive.mount(PATH)    # Google will require login credentials
    DATAPATH=PATH+'My Drive/data/'  # must end in "/"
    import requests
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_gen.py')
    with open('RNA_gen.py', 'w') as f:
        f.write(r.text)  # writes to cloud local, delete the file later?
    from RNA_gen import *
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_describe.py')
    with open('RNA_describe.py', 'w') as f:
        f.write(s.text)  # writes to cloud local, delete the file later?
    from RNA_describe import *
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_gen import *
    from SimTools.RNA_describe import *

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_gen():
    print("ERROR: Cannot use RNA_gen.")


# In[2]:


# Use code from our SimTools library.
def make_generator(seq_len):
    return cgen
def make_seqs(cgen,is_pc,train_count,test_count):
    freqs = [1,1,1,1]  # the relative frequencies for four nucleotides
    if is_pc:
        freqs = [2,1,1,2]  # protein-coding has more A and T
    else:
        pass # non-coding is random uniform
    cgen.get_seq_oracle().set_frequencies(freqs)    
    train_set = cgen.get_sequences(train_count)
    test_set =  cgen.get_sequences(test_count)
    return train_set,test_set


# In[3]:


import numpy as np
def get_the_facts(seqs):
    rd = RNA_describer()
    facts = rd.get_three_lengths(seqs)
    facts_ary = np.asarray(facts) # 5000 rows, 3 columns 
    print("Facts array:",type(facts_ary))
    print("Facts array:",facts_ary.shape)
    # Get the mean of each column
    mean_5utr, mean_orf, mean_3utr = np.mean(facts_ary,axis=0)
    std_5utr, std_orf, std_3utr = np.std(facts_ary,axis=0)
    print("mean 5' UTR length:",int(mean_5utr),"+/-",int(std_5utr))
    print("mean    ORF length:",int(mean_orf), "+/-",int(std_orf))
    print("mean 3' UTR length:",int(mean_3utr),"+/-",int(std_3utr))


# In[4]:


# Run this once to keep reusing same generator.
simulator = Collection_Generator() 
# If you want the same sequences every run, do this:
# simulator.set_reproducible(True) 


# ## Random uniform sequence

# In[9]:


# Run this repeatedly to get different random sequences.
SEQ_LEN=10000
SEQ_CNT=5000
simulator.get_len_oracle().set_mean(SEQ_LEN)
simulator.get_seq_oracle().set_sequences(['A','C','G','T'])  # default
simulator.get_seq_oracle().set_frequencies([1,1,1,1]) # default
seq_set = simulator.get_sequences(SEQ_CNT)
get_the_facts(seq_set)


# In[6]:


# In random sequence... 
#    any ORF tends to be short, so mean and std dev are small.
#    longest ORF tends toward 1st half of sequence (why?).
#    ORF can occur anywhere so UTR std dev is large.


# ## Non-random codon selection. 
# Expect longer ORFs than in purely random sequence.

# In[7]:


# Run this repeatedly to get different random sequences.
SEQ_LEN=333
SEQ_CNT=5000
simulator.get_len_oracle().set_mean(SEQ_LEN)
simulator.get_seq_oracle().set_sequences(['ATG','CCC','TAG'])
simulator.get_seq_oracle().set_frequencies([1,100,1]) 
seq_set = simulator.get_sequences(SEQ_CNT)
get_the_facts(seq_set)


# ## Simulated transcripts
# The simulator generates transcripts by the following formula.
# First 1/3 is random uniform letters to model 5'UTR.
# Middle 1/3 is an ORF: START codon, random non-stop codons, then STOP codon.
# Last 1/3 is random uniform letters to model 3'UTR.

# In[8]:


# Run this repeatedly to get different random sequences.
SEQ_LEN=1000
SEQ_CNT=5000
simulator.get_len_oracle().set_mean(SEQ_LEN)
simulator.set_seq_oracle(Transcript_Oracle())
seq_set = simulator.get_sequences(SEQ_CNT)
get_the_facts(seq_set)
# mean ORF length is slightly larger than 1/3 (why?)

