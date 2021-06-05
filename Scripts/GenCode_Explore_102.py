#!/usr/bin/env python
# coding: utf-8

# # GenCode Explore
# 
# Explore the human RNA sequences from GenCode.
# 
# Assume user downloaded files from GenCode 38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.
# 
# Exclude mitochondrial genes because many have non-standard start and stop codons.

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
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_special.py')
    with open('RNA_special.py', 'w') as f:
        f.write(s.text)  # writes to cloud local, delete the file later?
    from RNA_describe import *
    from RNA_special import *
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import *
    from SimTools.RNA_special import *

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# In[3]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'
TEST_FILENAME='test.fa.gz'


# In[4]:


def load_gencode(filename,label):
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
    special = RNA_Special_Cases()
    special.mitochondria()
    # Use gzip 'r' mode to open file in read-only mode.
    # Use gzip 't' mode to read each line of text as type string.
    with gzip.open (filename,'rt') as infile:
        for line in infile:
            if line[0]==DEFLINE:
                # Save the previous sequence (if previous exists).
                if not one_seq == EMPTY and not special.is_special(one_id):
                    labels.append(label)
                    seqs.append(one_seq)
                    lens.append(len(one_seq))
                    ids.append(one_id)
                # Get ready to read the next sequence. 
                # Parse a GenCode defline that is formatted like this:
                # >ENST0001234.5|gene_ID|other_fields other_info|other_info
                # Use the following if ever GenCode includes an ID with no version
                # one_id = line[1:].split(DELIM)[0].split(VERSION)[0] 
                one_id = line[1:].split(VERSION)[0]
                one_seq = EMPTY
            else:
                # Continue loading sequence lines till next defline.
                additional = line.rstrip()
                one_seq = one_seq + additional
        # Don't forget to save the last sequence after end-of-file.
        if not one_seq == EMPTY and not special.is_special(one_id):
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
pcdf=load_gencode(PC_FULLPATH,1)
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


def get_the_facts(seqs,verbose=False):
    oc = ORF_counter()
    count = len(seqs)
    max_orf_lengths=np.zeros(count)
    for s in range(0,count):
        seq = seqs[s]
        oc.set_sequence(seq)
        max_orf = oc.get_max_orf_len()
        max_orf_lengths[s] = max_orf
    mean_max_orf = np.mean(max_orf_lengths,axis=0)
    std_max_orf = np.std(max_orf_lengths,axis=0)
    if verbose:
        print("mean longest ORF length:",int(mean_max_orf),"+/-",int(std_max_orf))
    return mean_max_orf


# In[8]:


# Warning: each get_the_facts() can take up to 5 minutes.
# It is basically a 3-deep nested loop: for each seq, for each start, for each stop.
# Usually run this on subsets, not the whole data set.
def big_summary():
    show_time()
    print("Protein Coding set:")
    pc_means = get_the_facts( pcdf['sequence'].tolist() ,True)
    show_time()
    print("Non Coding set:")
    nc_means = get_the_facts( ncdf['sequence'].tolist() ,True)
    show_time()
big_summary() # this is slow


# GenCode38  
# ```
# 2021-05-28 16:22:23 EDT  
# Protein Coding set:  
# Facts array: (106143, 3)  
# mean 5' UTR length: 261 +/- 339  
# mean    ORF length: 1136 +/- 1556  
# mean 3' UTR length: 897 +/- 1385  
# 2021-05-28 16:26:34 EDT  
# Non Coding set:  
# Facts array: (48752, 3)  
# mean 5' UTR length: 511 +/- 1344  
# mean    ORF length: 211 +/- 135  
# mean 3' UTR length: 606 +/- 1100  
# 2021-05-28 16:27:00 EDT  
# ```

# ## Subset by RNA length and analyze ORF lengths
# 

# In[9]:


# This is a fast way to slice if you have length thresholds.
mask = (ncdf['sequence'].str.len() < 1000)
subset = ncdf.loc[mask]
discard = get_the_facts( subset['sequence'].tolist() ,True)


# In[11]:


def show_divisions(df,divisions,label):
    total=len(df)
    step=total//divisions
    for i in range(0,total,step):
        subset = df[i:i+step]
        first_len=subset.iloc[0]['seqlen']
        last_len=subset.iloc[-1]['seqlen']
        print("-- ",label,"RNA lengths",first_len,"to",last_len)
        discard = get_the_facts( subset['sequence'].tolist() ,True)
show_divisions(ncdf,10,"NC")
print()
show_divisions(pcdf,10,"PC")


# In[ ]:




