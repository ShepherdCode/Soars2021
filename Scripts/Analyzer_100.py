#!/usr/bin/env python
# coding: utf-8

# ##Analyzer
# 
# Statistically and visually compare mRNA and lncRNA sequences from GenCode.v38.
# 
# Assume user downloaded files from GenCode38 [FTP](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/)
# to a subdirectory called data.

# ##Import Dependencies
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import gzip
from scipy.stats import chisquare, kstest, spearmanr
import scipy.stats as ss
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
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/GenCodeTools.py')
    with open ('GenCodeTools.py', 'w') as f:
      f.write(s.text)
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/plot_generator.py')
    with open('plot_generator.py', 'w') as f:
      f.write(s.text)
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/KmerTools.py')
    with open('KmerTools.py', 'w') as f:
      f.write(s.text)  
    from KmerTools import KmerTools
    from RNA_describe import *
    from GenCodeTools import *
    from plot_generator import *
except:
    print("CoLab not working. On my PC, use relative paths.")
    IN_COLAB = False
    DATAPATH='../data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs
    from SimTools.RNA_describe import *
    from SimTools.GenCodeTools import *
    from SimTools.plot_generator import *
    from SimTools.KmerTools import KmerTools

MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login

if not assert_imported_RNA_describe():
    print("ERROR: Cannot use RNA_describe.")


# ##Load GenCode Data
# Loads GenCode.v38 data.
# 
# Filters out mRNA sequences based on UTR check.

# In[ ]:


PC_FILENAME='gencode.v38.pc_transcripts.fa.gz'
NC_FILENAME='gencode.v38.lncRNA_transcripts.fa.gz'
PC_FULLPATH=DATAPATH+PC_FILENAME
NC_FULLPATH=DATAPATH+NC_FILENAME
loader=GenCodeLoader()
loader.set_label(1)
loader.set_check_list(None) 
loader.set_check_utr(True)
pcdf=loader.load_file(PC_FULLPATH)
print("PC seqs loaded:",len(pcdf))
loader.set_label(0)
loader.set_check_list(None)
loader.set_check_utr(False)
ncdf=loader.load_file(NC_FULLPATH)
print("NC seqs loaded:",len(ncdf))


# ##Process Sequences
# Sampling, binning, length constraints, etc.
# 

# In[ ]:


APPLY_SUBSET = True             #Option to subset the data
MINIMUM_SEQUENCE_LENGTH = 0     #Minimum inclusive length to filter out sequences by
MAXIMUM_SEQUENCE_LENGTH = 4000  #Maximum inclusive length to filter out sequences by
SAMPLE_FRACTION = 0.5           #What fraction of the GenCode data set to take a sample of
REPRODUCABILITY_SEED = 314159   #Use to reproduce random sampling


# In[ ]:


if APPLY_SUBSET:
  pcdf = pcdf.sample(frac=SAMPLE_FRACTION, random_state=REPRODUCABILITY_SEED)
  ncdf = ncdf.sample(frac=SAMPLE_FRACTION, random_state=REPRODUCABILITY_SEED)

  print('PC sample size:', len(pcdf))
  print('NC sample size:', len(ncdf))


# In[ ]:


def subset_list_by_len_bounds(input_list, min_len, max_len):
  return list(filter(lambda x: len(x) > min_len and len(x) <= max_len, input_list))


# In[ ]:


pc_sequences = pcdf['sequence'].tolist()
nc_sequences = ncdf['sequence'].tolist()

if APPLY_SUBSET:
  pc_sequences = subset_list_by_len_bounds(pc_sequences, MINIMUM_SEQUENCE_LENGTH, MAXIMUM_SEQUENCE_LENGTH)
  nc_sequences = subset_list_by_len_bounds(nc_sequences, MINIMUM_SEQUENCE_LENGTH, MAXIMUM_SEQUENCE_LENGTH)

  print('PC seqs in length range','('+str(MINIMUM_SEQUENCE_LENGTH),'-',str(MAXIMUM_SEQUENCE_LENGTH)+'):', len(pc_sequences))
  print('NC seqs in length range','('+str(MINIMUM_SEQUENCE_LENGTH),'-',str(MAXIMUM_SEQUENCE_LENGTH)+'):', len(nc_sequences))

#Garbage collection
pcdf = None
ncdf = None


# ##Generate Statistics

# 
# Using KmerTools to get the K-mer counts upto 3.
# It returns the value in Dictionary form. (Key-Value Pair)

# In[ ]:


MAX_K = 3
tool = KmerTools()
pc_counts = tool.make_dict_upto_K(MAX_K)
for sample in pc_sequences:
    tool.update_count_one_K(pc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(pc_counts,MAX_K)
print("PC counts:\n",pc_counts)
pc_freqs = tool.count_to_frequency(pc_counts,MAX_K)
print ("Frequency:\n",pc_freqs)

nc_counts = tool.make_dict_upto_K(MAX_K)
for sample in nc_sequences:
    tool.update_count_one_K(nc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(nc_counts,MAX_K)
print("NC counts:\n",nc_counts)
nc_freqs = tool.count_to_frequency(nc_counts,MAX_K)
print ("Frequency:\n",nc_freqs)


# This function takes a dictionary as a parameter. It checks the length of the key to determine if it is 1-mer, 2-mer or 3 mer and assign the values respectively. 

# In[ ]:


def get_stats(dict):
  one_mer = []
  one_mer_key = []
  two_mer= []
  two_mer_key=[]
  three_mer= []
  three_mer_key=[]
  for sequence in enumerate(dict.items()):
    if(len(sequence[1][0])==1):
      one_mer.append(sequence[1][1])
      one_mer_key.append(sequence[1][0])
    if(len(sequence[1][0])==2):
      two_mer.append(sequence[1][1])
      two_mer_key.append(sequence[1][0])
    if(len(sequence[1][0])==3):
      three_mer.append(sequence[1][1])
      three_mer_key.append(sequence[1][0])
  return one_mer_key, one_mer, two_mer_key,two_mer, three_mer_key, three_mer


# In[51]:


#Gets the list  of 1-mer, 2-mer and 3-mer counts along with their key
#In the order of Key-Value. 
pc_stats = get_stats(pc_freqs)
nc_stats = get_stats(nc_freqs)

#1-mer
one_mer_pc = np.asarray(pc_stats[1])
one_mer_nc = np.asarray(nc_stats[1])
#2-mer
two_mer_pc = np.asarray(pc_stats[3])
two_mer_nc = np.asarray(nc_stats[3])
#3-mer
three_mer_pc = np.asarray(pc_stats[5])
three_mer_nc = np.asarray(nc_stats[5])

#Keys that can be used as labels.
one_mer_keys = pc_stats[0]
two_mer_keys = pc_stats[2]
three_mer_keys = pc_stats[4]

print(two_mer_pc)


# In[ ]:


print(one_mer_keys)
print(one_mer_pc)
print(one_mer_nc)


# In[ ]:


oc = ORF_counter()

pc_max_orf_len = np.empty(1, dtype=object)
nc_max_orf_len = np.empty(1, dtype=object)
pc_max_orf_len[0] = np.zeros(len(pc_sequences))
nc_max_orf_len[0] = np.zeros(len(nc_sequences))

for i in range(len(pc_sequences)):
  oc.set_sequence(pc_sequences[i])
  pc_max_orf_len[0][i] = oc.get_max_orf_len()
for i in range(len(nc_sequences)):
  oc.set_sequence(nc_sequences[i])
  nc_max_orf_len[0][i] = oc.get_max_orf_len()


# Using Correlation to see if the total length of series have any relation to the length of ORF. 

# In[ ]:


tools = ORF_RE() #Uses Regular Expression from Professor Miller's RNA_describe.

total_length_of_sequence=[]
length_of_ORF=[]
for i in range(len(pc_sequences)):
  total_length_of_sequence.append(len(pc_sequences[i]))
  length_of_ORF.append(tools.get_three_lengths(pc_sequences[i])[1])
r = np.corrcoef(np.asarray(total_length_of_sequence), np.asarray(length_of_ORF))
print(r)
plt.scatter(np.asarray(total_length_of_sequence), np.asarray(length_of_ORF))
plt.show()


# In[ ]:


total_length_of_sequence=[]
length_of_ORF=[]
for i in range(len(nc_sequences)):
  total_length_of_sequence.append(len(nc_sequences[i]))
  length_of_ORF.append(tools.get_three_lengths(nc_sequences[i])[1])
r = np.corrcoef(np.asarray(total_length_of_sequence), np.asarray(length_of_ORF))
print(r)
plt.scatter(np.asarray(total_length_of_sequence), np.asarray(length_of_ORF))
plt.show()


# KS-Similarity Test
# Not sure how it works but lets give it a shot. 

# In[ ]:


x = kstest(two_mer_pc, "norm")
print(x)


# ## Spearman Rank Correlation

# In[40]:


list_pc = random.sample(pc_sequences,100)
list_nc = random.sample(nc_sequences, 100)
type(pc_sequences[0])


# In[55]:


MAX_K = 6
tool = KmerTools()
kmer_pc = []
kmer_nc=[]
for sample in list_pc:
    pc_counts = tool.make_dict_upto_K(MAX_K)
    #tool.update_count_one_K(pc_counts,MAX_K,sample,True)
    #tool.harvest_counts_from_K(pc_counts,MAX_K)
    tool.update_count_one_K(pc_counts,MAX_K,sample)
    pc_freqs = tool.count_to_frequency(pc_counts,MAX_K)
    pc_vals = list(pc_freqs.values())
    kmer_pc.append(pc_vals)
#print(kmer_pc)



for sample in list_nc:
    nc_counts = tool.make_dict_upto_K(MAX_K)
    #tool.update_count_one_K(nc_counts,MAX_K,sample,True)
    #tool.harvest_counts_from_K(nc_counts,MAX_K)
    tool.update_count_one_K(nc_counts,MAX_K,sample)
    nc_freqs = tool.count_to_frequency(nc_counts,MAX_K)
    nc_vals = list(nc_freqs.values())
    kmer_nc.append(nc_vals)

count =0
while(count<1):
  #pc_stats = get_stats(kmer_pc[count])
  #nc_stats = get_stats(kmer_nc[count])
  #3-mer

  #three_mer_pc = np.asarray(pc_stats[5])
  #three_mer_nc = np.asarray(nc_stats[5])
  #print(three_mer_pc)
  #print(three_mer_nc)

  coef, p = spearmanr (kmer_pc[0],kmer_pc[1])
  coef, p = spearmanr (kmer_nc[0],kmer_nc[1])

  print(kmer_pc[0])
  print('Spearmans correlation coefficient: %.3f' % coef)
  alpha = 0.05
  #if p > alpha:
	  #print('Samples are uncorrelated (fail to reject H0) p=%.3f' % p)
  #else:
	  #print('Samples are correlated (reject H0) p=%.3f' % p)
  count+=1


# In[ ]:


#Gets the list  of 1-mer, 2-mer and 3-mer counts along with their key
#In the order of Key-Value. 
pc_stats = get_stats(pc_freqs)
nc_stats = get_stats(nc_freqs)


#3-mer
three_mer_pc = np.asarray(pc_stats[5])
three_mer_nc = np.asarray(nc_stats[5])

#Keys that can be used as labels.
one_mer_keys = pc_stats[0]
two_mer_keys = pc_stats[2]
three_mer_keys = pc_stats[4]


# In[ ]:





# ##Results

# 

# In[ ]:


data_set_names = ['mRNA', 'lncRNA']


# In[ ]:


pg = PlotGenerator()
pg.set_text_options(45, 'right', 0, 'center', 12)


# In[ ]:


pg.set_text('1-Mer Frequencies', 'Mers', 'Frequency', one_mer_keys, None)
pg.bar_plot([one_mer_pc, one_mer_nc], data_set_names)

pg.set_text('2-Mer Frequencies', 'Mers', 'Frequency', two_mer_keys, None)
pg.bar_plot([two_mer_pc, two_mer_nc], data_set_names)

pg.set_figure_options(width=16)
pg.set_text_options(90, 'right', 0, 'center', 7)
pg.set_text('3-Mer Frequencies', 'Mers', 'Frequency', three_mer_keys, None)
pg.bar_plot([three_mer_pc, three_mer_nc], data_set_names)

pg.set_figure_options()
pg.set_text_options(45, 'right', 0, 'center', 12)
pg.set_text('Max ORF Lengths', 'RNA Types', 'Max ORF Length', [''], None)
pg.box_plot([pc_max_orf_len, nc_max_orf_len], data_set_names, False)

