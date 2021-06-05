#!/usr/bin/env python
# coding: utf-8

# # GenCode Explore
# 
# Explore the human RNA sequences from GenCode.
# 
# Assume user downloaded files from GenCode 38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.
# 
# Improve on GenCode_Explore_101.ipynb
# 
# Use ORF_counter. 
# 
# Use MatPlotLib to make box plots and heat maps.

# In[150]:


import time 
def show_time():
    t = time.time()
    s = time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t))
    print(s)
show_time()


# In[151]:


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
    from RNA_describe import *
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import *

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# In[152]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'


# In[153]:


def load_gencode(filename,label):
    DEFLINE='>'
    DELIM='|'
    EMPTY=''
    labels=[]  # usually 1 for protein-coding or 0 for non-coding
    seqs=[]    # usually string of ACGT
    lens=[]    # sequence length
    ids=[]     # GenCode transcript ID, always starts ENST
    one_seq = EMPTY
    one_id = None
    # Use gzip 'r' mode to open file in read-only mode.
    # Use gzip 't' mode to read each line of text as type string.
    with gzip.open (filename,'rt') as infile:
        for line in infile:
            if line[0]==DEFLINE:
                # Save the previous sequence if one exists.
                if not one_seq == EMPTY:
                    labels.append(label)
                    seqs.append(one_seq)
                    lens.append(len(one_seq))
                    ids.append(one_id)
                # Get ready to read the next sequence. 
                # Parse a GenCode defline that is formatted like this:
                # >transcript_ID|gene_ID|other_fields other_info|other_info
                one_id = line[1:].split(DELIM)[0]
                one_seq = EMPTY
            else:
                # Continue loading sequence lines till next defline.
                additional = line.rstrip()
                one_seq = one_seq + additional
        # Don't forget to save the last sequence after end-of-file.
        if not one_seq == EMPTY:
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


# In[154]:


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


# ## Load the GenCode data.
# Warning: GenCode has
# over 100K protein-coding RNA (mRNA) 
# and almost 50K non-coding RNA (lncRNA).

# In[155]:


PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME
show_time()
pcdf=load_gencode(PC_FULLPATH,1)
print("PC seqs loaded:",len(pcdf))
show_time()
ncdf=load_gencode(NC_FULLPATH,0)
print("NC seqs loaded:",len(ncdf))
show_time()


# In[156]:


print("Sorting PC...")
pcdf.sort_values('seqlen', ascending=True, inplace=True)
print("Sorting NC...")
ncdf.sort_values('seqlen', ascending=True, inplace=True)


# ## Subset by RNA length and analyze ORF lengths
# 
# ---
# 
# 
# 

# In[157]:


# This is a fast way to slice if you have length thresholds.
# TO DO: choose length thresholds and apply to PC and NC RNA.
# For example: 200, 400, 800, 1600, 3200, 6400 (e.g. 200-399, etc.)
mask = (ncdf['sequence'].str.len() < 1000)
subset = ncdf.loc[mask]

# Here is one way to extract a list from a dataframe. 
mylist=subset['sequence'].tolist()

mask = (pcdf['sequence'].str.len() < 800)
subset = pcdf.loc[mask]
subset_list = subset['sequence'].tolist()


# In[158]:


def subset_list_by_len_bounds(input_list, min_len, max_len):
  return list(filter(lambda x: len(x) > min_len and len(x) < max_len, input_list))


# In[159]:


import matplotlib.pyplot as plt
import numpy as np

#TODO: also do non-coding sequences
#TODO: optimize the application of length thresholds and make use of numpy lists

x_axis_labels = []

pc_max_len_data = []
pc_max_cnt_data = []
pc_contain_data = []

for bin_num in range(1, 7):
  pc_max_len_data.append([])
  pc_max_cnt_data.append([])
  pc_contain_data.append([])

  #define min and max sequence ranges (200, 400, 800, 1600, etc.)
  min_seq_length = 2**bin_num * 100
  max_seq_length = 2**(bin_num + 1) * 100

  #just for the box plots
  x_axis_labels.append(str(min_seq_length) + "-" + str(max_seq_length))

  #subset the protein-coding sequences by min and max sequence length bounds
  pc_sequences = pcdf['sequence'].tolist()
  pc_sequences_subset = subset_list_by_len_bounds(pc_sequences, min_seq_length, max_seq_length)

  #get stats on each protein-coding sequence
  for seq_num in range(0, len(pc_sequences_subset)):
    sequence = pc_sequences_subset[seq_num]
    oc = ORF_counter()
    oc.set_sequence(sequence)

    max_len = oc.get_max_orf_len()
    max_cnt = oc.count_maximal_orfs()
    contain = oc.count_contained_orfs()

    pc_max_len_data[bin_num - 1].append(max_len)
    pc_max_cnt_data[bin_num - 1].append(max_cnt)
    pc_contain_data[bin_num - 1].append(contain)


# ## Plot the Data
# 
# ---

# In[160]:


#plot stats of protein-coding sequences
fig1, ax1 = plt.subplots()
ax1.set_yscale('log', basey=2)
ax1.set_title('Length of Longest ORF in Protein-coding Sequences')
ax1.set_ylabel('ORF Length')
ax1.set_xlabel('Sequence Length Ranges')
ax1.set_xticklabels(x_axis_labels, rotation = 45, ha = "center")
ax1.boxplot(pc_max_len_data, showfliers=False)

fig1, ax1 = plt.subplots()
ax1.set_yscale('log', basey=2)
ax1.set_title('Number of Non-contained ORFs in Protein-coding Sequences')
ax1.set_ylabel('Number of Non-contained ORFs')
ax1.set_xlabel('Sequence Length Ranges')
ax1.set_xticklabels(x_axis_labels, rotation = 45, ha = "center")
ax1.boxplot(pc_max_cnt_data, showfliers=False)

fig1, ax1 = plt.subplots()
ax1.set_yscale('log', basey=2)
ax1.set_title('Number of Contained ORFs in Protein-coding Sequences')
ax1.set_ylabel('Number of Contained ORFs')
ax1.set_xlabel('Sequence Length Ranges')
ax1.set_xticklabels(x_axis_labels, rotation = 45, ha = "center")
ax1.boxplot(pc_contain_data, showfliers=False)


# ## Plotting examples
# [boxplot doc](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html)  
# [boxplot demo](https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html)  
# [heatmap examples](https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib) - scroll down!  

# In[160]:




