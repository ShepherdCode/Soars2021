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

# In[38]:


import time 
def show_time():
    t = time.time()
    s = time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t))
    print(s)
show_time()


# In[39]:


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


# In[40]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'


# In[41]:


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


# In[42]:


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

# In[43]:


PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME
show_time()
pcdf=load_gencode(PC_FULLPATH,1)
print("PC seqs loaded:",len(pcdf))
show_time()
ncdf=load_gencode(NC_FULLPATH,0)
print("NC seqs loaded:",len(ncdf))
show_time()


# In[44]:


print("Sorting PC...")
pcdf.sort_values('seqlen', ascending=True, inplace=True)
print("Sorting NC...")
ncdf.sort_values('seqlen', ascending=True, inplace=True)


# ## Subset by RNA length
# 
# ---
# 
# 
# 

# In[45]:


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


# In[46]:


def subset_list_by_len_bounds(input_list, min_len, max_len):
  return list(filter(lambda x: len(x) > min_len and len(x) < max_len, input_list))


# In[47]:


import matplotlib.pyplot as plt
import numpy as np

x_axis_labels = []

#Bin the RNA sequences
bins = [(200, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 6400), (6400, 12800)]
binned_pc_sequences = []
binned_nc_sequences = []
for i in range(0, len(bins)):
  bin = bins[i]
  binned_pc_sequences.append([])
  binned_nc_sequences.append([])
  binned_pc_sequences[i] = subset_list_by_len_bounds(pcdf['sequence'].tolist(), bin[0], bin[1])
  binned_nc_sequences[i] = subset_list_by_len_bounds(ncdf['sequence'].tolist(), bin[0], bin[1])

  #Labels for the plots
  x_axis_labels.append(str(bin[0]) + "-" + str(bin[1]) + " (mRNA)")
  x_axis_labels.append(str(bin[0]) + "-" + str(bin[1]) + " (lncRNA)")


# ##Gather data on ORF lengths and the number of contained and non-contained ORFs
# 
# ---

# In[48]:


pc_max_len_data = []
pc_max_cnt_data = []
pc_contain_data = []

nc_max_len_data = []
nc_max_cnt_data = []
nc_contain_data = []

oc = ORF_counter()

for bin_num in range(0, len(bins)):

  pc_max_len_data.append([])
  pc_max_cnt_data.append([])
  pc_contain_data.append([])
  for seq_num in range(0, len(binned_pc_sequences[bin_num])):
    oc.set_sequence(binned_pc_sequences[bin_num][seq_num])
    pc_max_len_data[bin_num].append(oc.get_max_orf_len())
    pc_max_cnt_data[bin_num].append(oc.count_maximal_orfs())
    pc_contain_data[bin_num].append(oc.count_contained_orfs())

  nc_max_len_data.append([])
  nc_max_cnt_data.append([])
  nc_contain_data.append([])
  for seq_num in range(0, len(binned_nc_sequences[bin_num])):
    oc.set_sequence(binned_nc_sequences[bin_num][seq_num])
    nc_max_len_data[bin_num].append(oc.get_max_orf_len())
    nc_max_cnt_data[bin_num].append(oc.count_maximal_orfs())
    nc_contain_data[bin_num].append(oc.count_contained_orfs())


# ##Merge protein-coding and non-coding data for plotting
# 
# ---

# In[49]:


max_len_data = []
max_cnt_data = []
contain_data = []
for i in range(0, len(bins)):
  max_len_data.append(pc_max_len_data[i])
  max_len_data.append(nc_max_len_data[i])
  max_cnt_data.append(pc_max_cnt_data[i])
  max_cnt_data.append(nc_max_cnt_data[i])
  contain_data.append(pc_contain_data[i])
  contain_data.append(nc_contain_data[i])


# ## Plot the data
# 
# ---

# In[50]:


fig1, ax1 = plt.subplots()
ax1.set_yscale('log', basey=2)
ax1.set_title('Length of Longest ORF in RNA Sequences')
ax1.set_ylabel('ORF Length')
ax1.set_xlabel('Sequence Length Ranges')
ax1.set_xticklabels(x_axis_labels, rotation = 45, ha = "right")
ax1.boxplot(max_len_data, showfliers=False)

fig2, ax2 = plt.subplots()
ax2.set_yscale('log', basey=2)
ax2.set_title('Number of Non-contained ORFs in RNA Sequences')
ax2.set_ylabel('Number of Non-contained ORFs')
ax2.set_xlabel('Sequence Length Ranges')
ax2.set_xticklabels(x_axis_labels, rotation = 45, ha = "right")
ax2.boxplot(max_cnt_data, showfliers=False)

fig3, ax3 = plt.subplots()
ax3.set_yscale('log', basey=2)
ax3.set_title('Number of Contained ORFs in RNA Sequences')
ax3.set_ylabel('Number of Contained ORFs')
ax3.set_xlabel('Sequence Length Ranges')
ax3.set_xticklabels(x_axis_labels, rotation = 45, ha = "right")
ax3.boxplot(contain_data, showfliers=False)


# ## Plotting examples
# [boxplot doc](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html)  
# [boxplot demo](https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html)  
# [heatmap examples](https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib) - scroll down!  

# In[50]:




