#!/usr/bin/env python
# coding: utf-8

# # GenCode Explore
# 
# Explore the human RNA sequences from GenCode.
# 
# Assume user downloaded files from GenCode 38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.
# 
# Move the GenCodeLoader class to its own python module. Compare to 105.

# In[1]:


import time 
def show_time():
    t = time.time()
    s = time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t))
    print(s)
show_time()


# In[2]:


import numpy as np
import pandas as pd
import sys

try:
    from google.colab import drive
    IN_COLAB = True
    print("On Google CoLab, mount cloud-local file, get our code from GitHub.")
    PATH='/content/drive/'
    #drive.mount(PATH,force_remount=True)  # hardly ever need this
    drive.mount(PATH)    # Google will require login credentials
    DATAPATH=PATH+'My Drive/data/'  # must end in "/"
    import requests
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_describe.py')
    with open('RNA_describe.py', 'w') as f:
        f.write(s.text)  # writes to cloud local, delete the file later?
    from RNA_describe import ORF_counter
    from RNA_describe import assert_imported_RNA_describe
    from GenCodeTools import GenCodeLoader
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import ORF_counter
    from SimTools.RNA_describe import assert_imported_RNA_describe
    from SimTools.GenCodeTools import GenCodeLoader

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# In[3]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'


# ## Load the GenCode data.
# Warning: GenCode has
# over 100K protein-coding RNA (mRNA) 
# and almost 50K non-coding RNA (lncRNA).

# In[4]:


# Full GenCode ver 38 human is 106143 pc + 48752 nc and loads in 7 sec.
# Expect fewer transcripts if special filtering is used.
PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME
loader=GenCodeLoader()
show_time()
loader.set_label(1)
loader.set_check_list(None) 
loader.set_check_utr(True)
pcdf=loader.load_file(PC_FULLPATH)
print("PC seqs loaded:",len(pcdf))
show_time()
loader.set_label(0)
loader.set_check_list(None)
loader.set_check_utr(False)
ncdf=loader.load_file(NC_FULLPATH)
print("NC seqs loaded:",len(ncdf))
show_time()


# In[5]:


print("Sorting PC...")
pcdf.sort_values('seqlen', ascending=True, inplace=True)
print("Sorting NC...")
ncdf.sort_values('seqlen', ascending=True, inplace=True)


# In[6]:


ncdf


# ## Look for short ORFs

# In[7]:


def show_short(df,too_short):
    oc = ORF_counter()
    count=len(df)
    shorties=0
    for pos in range(0,count):
        sequence=df.iloc[pos]['sequence']
        seqlen=df.iloc[pos]['seqlen']
        oc.set_sequence(sequence)
        orflen=oc.get_max_orf_len()
        seqlen=df.iloc[pos]['seqlen']
        if seqlen>200 and orflen<=TOO_SHORT:
            seqid=df.iloc[pos]['tid']
            #print("%s len=%d orf=%d"%(seqid,seqlen,orflen))
            shorties += 1
        if pos%10000==0:
            print("Up to position %d, we have %d shorter than %d"%(pos,shorties,too_short))
    print("After all %d, we have %d shorter than %d"%(count,shorties,too_short))
TOO_SHORT=60
show_short(pcdf,TOO_SHORT)


# In[8]:


show_short(ncdf,TOO_SHORT)


# ## Conclusion
# With TOO_SHORT=30 
# NON-CODING
# We have 589 shorter than 30, with most of them (504) shorter than 10000
# 
# CODING
# Using check_utr and check_list on pcdf, we have 0 shorter than 30.
# Using check_utr only, we have 0 shorter than 30.
# 
