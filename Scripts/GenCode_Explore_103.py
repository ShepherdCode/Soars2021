#!/usr/bin/env python
# coding: utf-8

# # GenCode Explore
# 
# Explore the human RNA sequences from GenCode.
# 
# Assume user downloaded files from GenCode 38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.
# 
# Build on 102 which excluded mitochondrial genes by ID. (If the number of exclusions grow, may need to move from an exclusion list to an annotation gff parser.)
# 
# Explore remaining PC mRNA that have tiny ORFs.

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
import gzip
import sys
import re

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
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/GenCode_Protein_Include.py')
    with open('GenCode_Protein_Include', 'w') as f:
        f.write(s.text)  # writes to cloud local, delete the file later?
    from RNA_describe import *
    from GenCode_preprocess import prot_incl
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import *
    from SimTools.GenCode_Protein_Include import prot_incl

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# In[3]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'


# In[4]:


def load_gencode(filename,label,check_list=None):
    DEFLINE='>'  # start of line with ids in a FASTA FILE
    DELIM='|'    # character between ids
    VERSION='.'  # character between id and version
    EMPTY=''     # use this to avoid saving "previous" sequence in first iteration
    labels=[]  # usually 1 for protein-coding or 0 for non-coding
    seqs=[]    # usually strings of ACGT
    lens=[]    # sequence length
    ids=[]     # GenCode transcript ID, always starts ENST, excludes version
    one_seq = EMPTY
    one_id = None
    pattern5=re.compile('.*UTR5:')
    pattern3=re.compile('.*UTR3:')
    has_utr=False
    with gzip.open (filename,'rt') as infile:
        for line in infile:
            if line[0]==DEFLINE:
                if not one_seq == EMPTY and                 (check_list is None or one_id in check_list) and has_utr:
                    labels.append(label)
                    seqs.append(one_seq)
                    lens.append(len(one_seq))
                    ids.append(one_id)
                one_id = line[1:].split(VERSION)[0]
                one_seq = EMPTY
                has_utr = not (pattern5.match(line) is None or pattern3.match(line) is None)
            else:
                # Continue loading sequence lines till next defline.
                additional = line.rstrip()
                one_seq = one_seq + additional
        # Don't forget to save the last sequence after end-of-file.
        if not one_seq == EMPTY and (check_list is None or one_id in check_list):
            labels.append(label)
            seqs.append(one_seq)
            lens.append(len(one_seq))
            ids.append(one_id)

    df1=pd.DataFrame(ids,columns=['tid'])
    df2=pd.DataFrame(labels,columns=['class'])
    df3=pd.DataFrame(seqs,columns=['sequence'])
    df4=pd.DataFrame(lens,columns=['seqlen'])
    df=pd.concat((df1,df2,df3,df4),axis=1)
    return df


# ## Load the GenCode data.
# Warning: GenCode has
# over 100K protein-coding RNA (mRNA) 
# and almost 50K non-coding RNA (lncRNA).

# In[5]:


# Full GenCode ver 38 human is 106143 pc + 48752 nc and loads in 7 sec.
# Expect fewer transcripts if special filtering is used.
PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME
show_time()
pcdf=load_gencode(PC_FULLPATH,1,prot_incl)
print("PC seqs loaded:",len(pcdf))
show_time()
ncdf=load_gencode(NC_FULLPATH,0)
print("NC seqs loaded:",len(ncdf))
show_time()


# In[6]:


print("Sorting PC...")
pcdf.sort_values('seqlen', ascending=True, inplace=True)
print("Sorting NC...")
ncdf.sort_values('seqlen', ascending=True, inplace=True)


# In[7]:


ncdf


# ## Look for short ORFs

# In[10]:


def show_short(df,too_short):
    oc = ORF_counter()
    count=len(df)
    for pos in range(0,count):
        sequence=df.iloc[pos]['sequence']
        seqlen=df.iloc[pos]['seqlen']
        oc.set_sequence(sequence)
        orflen=oc.get_max_orf_len()
        seqlen=df.iloc[pos]['seqlen']
        if seqlen>200 and orflen<=TOO_SHORT:
            seqid=df.iloc[pos]['tid']
            print("%s len=%d orf=%d"%(seqid,seqlen,orflen))
        if pos%10000==0:
            print("...up to position",pos)
    print("done")
TOO_SHORT=50
show_short(pcdf,TOO_SHORT)


# In[ ]:





# In[ ]:




