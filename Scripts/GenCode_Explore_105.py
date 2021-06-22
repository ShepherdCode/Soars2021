#!/usr/bin/env python
# coding: utf-8

# # GenCode Explore
# 
# Explore the human RNA sequences from GenCode.
# 
# Assume user downloaded files from GenCode 38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.
# 
# In 104, we showed that we can do away with the protein-include file based on annotation.gff and just rely on the presence of UTR in the FASTA deflines. Here, stop importing the protein-include file.

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
    from RNA_describe import ORF_counter
    from RNA_describe import assert_imported_RNA_describe
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import ORF_counter
    from SimTools.RNA_describe import assert_imported_RNA_describe

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# In[3]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'


# In[4]:


class GenCodeLoader():
    def __init__(self):
        self.pattern5=re.compile('.*UTR5:')
        self.pattern3=re.compile('.*UTR3:')
        self.check_list = None
        self.check_utr = False
    def set_label(self,label):
        self.label=label
    def set_check_list(self,check_list):
        self.check_list=check_list
    def set_check_utr(self,check_utr):
        self.check_utr=check_utr
    def __save_previous(self,one_def,one_seq):
        if one_def is None:
            return
        if self.check_utr:
            if self.pattern5.match(one_def) is None: 
                return
            if self.pattern3.match(one_def) is None:
                return
        VERSION = '.'
        one_id = one_def[1:].split(VERSION)[0]
        if self.check_list is not None:
            if one_id not in self.check_list:
                return
        self.labels.append(self.label)
        self.seqs.append(one_seq)
        self.lens.append(len(one_seq))
        self.ids.append(one_id)
    def load_file(self,filename):
        self.labels=[]  # usually 1 for protein-coding or 0 for non-coding
        self.seqs=[]    # usually strings of ACGT
        self.lens=[]    # sequence length
        self.ids=[]     # GenCode transcript ID, always starts ENST, excludes version
        DEFLINE='>'  # start of line with ids in a FASTA FILE
        EMPTY=''
        one_def = None
        one_seq = ''
        with gzip.open (filename,'rt') as infile:
            for line in infile:
                if line[0]==DEFLINE:
                    self.__save_previous(one_def,one_seq)
                    one_def=line
                    one_seq = EMPTY
                else:
                    # Continue loading sequence lines till next defline.
                    additional = line.rstrip()
                    one_seq = one_seq + additional
            # Don't forget to save the last sequence after end-of-file.
            self.__save_previous(one_def,one_seq)
        df1=pd.DataFrame(self.ids,columns=['tid'])
        df2=pd.DataFrame(self.labels,columns=['class'])
        df3=pd.DataFrame(self.seqs,columns=['sequence'])
        df4=pd.DataFrame(self.lens,columns=['seqlen'])
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


# In[6]:


print("Sorting PC...")
pcdf.sort_values('seqlen', ascending=True, inplace=True)
print("Sorting NC...")
ncdf.sort_values('seqlen', ascending=True, inplace=True)


# In[7]:


ncdf


# ## Look for short ORFs

# In[8]:


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


# In[9]:


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
