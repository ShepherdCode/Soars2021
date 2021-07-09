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
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/GenCodeTools.py')
    with open ('GenCodeTools.py', 'w') as f:
      f.write(s.text)
    s = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/plot_generator.py')
    with open('plot_generator.py', 'w') as f:
      f.write(s.text)
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


# ###OPTIONS
# ---

# In[5]:


SAMPLE_FRACTION = 0.5
REPRODUCABILITY_SEED = 314159 #Not used in Random_Base_Oracle or PlotGenerator
NUM_BINS = 3
BINNING_SCALE = 100
NUMBER_OF_SIM_SEQUENCES_PER_BIN = 100


# ###Take sample of mRNA and lncRNA sequence data sets
# ---

# In[6]:


pcdf_sample = pcdf.sample(frac=SAMPLE_FRACTION, random_state=REPRODUCABILITY_SEED) #Take a sample of a fraction of the mRNA data frame
ncdf_sample = ncdf.sample(frac=SAMPLE_FRACTION, random_state=REPRODUCABILITY_SEED) #Take a sample of a fraction of the lncRNA data frame
print('Sample size of mRNA sequences:', len(pcdf_sample))
print('Sample size of lncRNA sequences:', len(ncdf_sample))
show_time()


# ###Generate bins
# ---

# In[7]:


bins = []
for b in range(1, NUM_BINS + 1):
  bin = (2 ** (b) * BINNING_SCALE, 2 ** (b + 1) * BINNING_SCALE)
  bins.append(bin)
print('Bins generated:', bins)
show_time()


# ###Generate random RNA sequences
# 
# ---

# In[8]:


random.seed(REPRODUCABILITY_SEED)
sim_pc_sequences = []
sim_nc_sequences = []
for bin in bins:
  seq_len = (bin[0] + bin[1]) >> 1 #integer division using // apparently still calculates remainder, source: https://stackoverflow.com/questions/37053379/times-two-faster-than-bit-shift-for-python-3-x-integers
  cds_len = seq_len >> 1 #So instead use right bit shift
  print('Sequence length:',seq_len)
  rbo = Random_Base_Oracle(seq_len, debug=True)
  temp_pc, temp_nc = rbo.get_partitioned_sequences(cds_len, NUMBER_OF_SIM_SEQUENCES_PER_BIN)
  sim_pc_sequences += temp_pc
  sim_nc_sequences += temp_nc
print('mRNA sequences generated:', NUMBER_OF_SIM_SEQUENCES_PER_BIN * NUM_BINS)
print('lncRNA sequences generated:', NUMBER_OF_SIM_SEQUENCES_PER_BIN * NUM_BINS)
show_time()


# ###Bin sequences by sequence length
# ---

# In[9]:


def subset_list_by_len_bounds(input_list, min_len, max_len):
  return list(filter(lambda x: len(x) > min_len and len(x) <= max_len, input_list))


# In[10]:


#Bin the RNA sequences
binned_pc_sequences = []
binned_nc_sequences = []
binned_sim_pc_sequences = []
binned_sim_nc_sequences = []
for i in range(0, NUM_BINS):
  bin = bins[i]
  binned_pc_sequences.append(subset_list_by_len_bounds(pcdf_sample['sequence'].tolist(), bin[0], bin[1]))
  binned_nc_sequences.append(subset_list_by_len_bounds(ncdf_sample['sequence'].tolist(), bin[0], bin[1]))
  binned_sim_pc_sequences.append(subset_list_by_len_bounds(sim_pc_sequences, bin[0], bin[1]))
  binned_sim_nc_sequences.append(subset_list_by_len_bounds(sim_nc_sequences, bin[0], bin[1]))
show_time()


# ##Gather data on ORF lengths and the number of contained and non-contained ORFs
# 
# ---

# In[11]:


#TODO: optimize. combine data?
pc_max_len_data = np.empty(NUM_BINS, dtype=object)
pc_max_cnt_data = np.empty(NUM_BINS, dtype=object)
pc_contain_data = np.empty(NUM_BINS, dtype=object)
nc_max_len_data = np.empty(NUM_BINS, dtype=object)
nc_max_cnt_data = np.empty(NUM_BINS, dtype=object)
nc_contain_data = np.empty(NUM_BINS, dtype=object)
sim_pc_max_len_data = np.empty(NUM_BINS, dtype=object)
sim_pc_max_cnt_data = np.empty(NUM_BINS, dtype=object)
sim_pc_contain_data = np.empty(NUM_BINS, dtype=object)
sim_nc_max_len_data = np.empty(NUM_BINS, dtype=object)
sim_nc_max_cnt_data = np.empty(NUM_BINS, dtype=object)
sim_nc_contain_data = np.empty(NUM_BINS, dtype=object)

oc = ORF_counter()

for bin in range(0, NUM_BINS):
  pc_max_len_data[bin] = np.zeros(len(binned_pc_sequences[bin]))
  pc_max_cnt_data[bin] = np.zeros(len(binned_pc_sequences[bin]))
  pc_contain_data[bin] = np.zeros(len(binned_pc_sequences[bin]))
  nc_max_len_data[bin] = np.zeros(len(binned_nc_sequences[bin]))
  nc_max_cnt_data[bin] = np.zeros(len(binned_nc_sequences[bin]))
  nc_contain_data[bin] = np.zeros(len(binned_nc_sequences[bin]))
  sim_pc_max_len_data[bin] = np.zeros(len(binned_sim_pc_sequences[bin]))
  sim_pc_max_cnt_data[bin] = np.zeros(len(binned_sim_pc_sequences[bin]))
  sim_pc_contain_data[bin] = np.zeros(len(binned_sim_pc_sequences[bin]))
  sim_nc_max_len_data[bin] = np.zeros(len(binned_sim_nc_sequences[bin]))
  sim_nc_max_cnt_data[bin] = np.zeros(len(binned_sim_nc_sequences[bin]))
  sim_nc_contain_data[bin] = np.zeros(len(binned_sim_nc_sequences[bin]))
  #Gather protein-coding sequence data
  for seq in range(0, len(binned_pc_sequences[bin])):
    oc.set_sequence(binned_pc_sequences[bin][seq])
    pc_max_len_data[bin][seq] = oc.get_max_orf_len()
    pc_max_cnt_data[bin][seq] = oc.count_maximal_orfs()
    pc_contain_data[bin][seq] = oc.count_contained_orfs()

  #Gather non-coding sequence data
  for seq in range(0, len(binned_nc_sequences[bin])):
    oc.set_sequence(binned_nc_sequences[bin][seq])
    nc_max_len_data[bin][seq] = oc.get_max_orf_len()
    nc_max_cnt_data[bin][seq] = oc.count_maximal_orfs()
    nc_contain_data[bin][seq] = oc.count_contained_orfs()

  #Gather simulated pc sequence data
  for seq in range(0, len(binned_sim_pc_sequences[bin])):
    oc.set_sequence(binned_sim_pc_sequences[bin][seq])
    sim_pc_max_len_data[bin][seq] = oc.get_max_orf_len()
    sim_pc_max_cnt_data[bin][seq] = oc.count_maximal_orfs()
    sim_pc_contain_data[bin][seq] = oc.count_contained_orfs()

  #Gather simulated nc sequence data
  for seq in range(0, len(binned_sim_nc_sequences[bin])):
    oc.set_sequence(binned_sim_nc_sequences[bin][seq])
    sim_nc_max_len_data[bin][seq] = oc.get_max_orf_len()
    sim_nc_max_cnt_data[bin][seq] = oc.count_maximal_orfs()
    sim_nc_contain_data[bin][seq] = oc.count_contained_orfs()
show_time()


# ###Calculate Mean Base Frequencies
# 
# ---

# In[14]:


BASES = ['A', 'C', 'G', 'T']
def count_bases(sequence):
  """
  Count the bases in a sequence.
  Returns a list where [0]=A, [1]=C, [2]=G, and [3]=T.
  """
  base_counts = np.zeros(4)
  for base in sequence:
    for i in range(0, len(BASES)):
      if base == BASES[i]:
        base_counts[i] += 1
        break
  return base_counts

def base_frequencies(sequence):
  """
  Calculate the frequency of each base in a sequence.
  Returns a list where [0]=A, [1]=C, [2]=G, and [3]=T.
  """
  return count_bases(sequence) / len(sequence)


# In[15]:


mean_pc_base_frequencies = np.zeros(4)
mean_nc_base_frequencies = np.zeros(4)
mean_sim_pc_base_frequencies = np.zeros(4)
mean_sim_nc_base_frequencies = np.zeros(4)
num_pc_seqs = 0
num_nc_seqs = 0
num_sim_pc_seqs = 0
num_sim_nc_seqs = 0
for i in range(0, NUM_BINS):
  for seq in binned_pc_sequences[i]:
    mean_pc_base_frequencies += base_frequencies(seq)
    num_pc_seqs += 1
  for seq in binned_nc_sequences[i]:
    mean_nc_base_frequencies += base_frequencies(seq)
    num_nc_seqs += 1
  for seq in binned_sim_pc_sequences[i]:
    mean_sim_pc_base_frequencies += base_frequencies(seq)
    num_sim_pc_seqs += 1
  for seq in binned_sim_nc_sequences[i]:
    mean_sim_nc_base_frequencies += base_frequencies(seq)
    num_sim_nc_seqs += 1
mean_pc_base_frequencies = mean_pc_base_frequencies / num_pc_seqs
mean_nc_base_frequencies = mean_nc_base_frequencies / num_nc_seqs
mean_sim_pc_base_frequencies = mean_sim_pc_base_frequencies / num_sim_pc_seqs
mean_sim_nc_base_frequencies = mean_sim_nc_base_frequencies / num_sim_nc_seqs


# ###Prepare data for plot of bin sizes
# 
# ---

# In[16]:


pc_bin_sizes = np.zeros(NUM_BINS)
nc_bin_sizes = np.zeros(NUM_BINS)
sim_pc_bin_sizes = np.zeros(NUM_BINS)
sim_nc_bin_sizes = np.zeros(NUM_BINS)
for i in range(0, NUM_BINS):
  pc_bin_sizes[i] = len(binned_pc_sequences[i])
  nc_bin_sizes[i] = len(binned_nc_sequences[i])
  sim_pc_bin_sizes[i] = len(binned_sim_pc_sequences[i])
  sim_nc_bin_sizes[i] = len(binned_sim_nc_sequences[i])
show_time()


# ###Prepare data for plot of number of sequences with no ORFs and plot of number of sequences with max ORF lengths equal to or less than 100
# 
# ---

# In[17]:


"""
Count the number of values in a given data set that are within a given inclusive range.
"""
def count_data_in_range(data, min, max):
  return np.sum((data >= min) & (data <= max))


# In[18]:


pc_no_orf_count = np.zeros(NUM_BINS)
nc_no_orf_count = np.zeros(NUM_BINS)
sim_pc_no_orf_count = np.zeros(NUM_BINS)
sim_nc_no_orf_count = np.zeros(NUM_BINS)
pc_max_orf_len_less_than_100 = np.zeros(NUM_BINS)
nc_max_orf_len_less_than_100 = np.zeros(NUM_BINS)
sim_pc_max_orf_len_less_than_100 = np.zeros(NUM_BINS)
sim_nc_max_orf_len_less_than_100 = np.zeros(NUM_BINS)
for i in range(0, NUM_BINS):
  pc_no_orf_count[i] = count_data_in_range(pc_max_len_data[i], 0, 0)
  nc_no_orf_count[i] = count_data_in_range(nc_max_len_data[i], 0, 0)
  sim_pc_no_orf_count[i] = count_data_in_range(sim_pc_max_len_data[i], 0, 0)
  sim_nc_no_orf_count[i] = count_data_in_range(sim_nc_max_len_data[i], 0, 0)
  pc_max_orf_len_less_than_100[i] = count_data_in_range(pc_max_len_data[i], 0, 100)
  nc_max_orf_len_less_than_100[i] = count_data_in_range(nc_max_len_data[i], 0, 100)
  sim_pc_max_orf_len_less_than_100[i] = count_data_in_range(sim_pc_max_len_data[i], 0, 100)
  sim_nc_max_orf_len_less_than_100[i] = count_data_in_range(sim_nc_max_len_data[i], 0, 100)

show_time()


# ## Plot the data
# 
# ---

# In[19]:


#Generate x-axis labels
x_axis_labels = []
for bin in bins:
  x_axis_labels.append(str(bin[0]) + "-" + str(bin[1]))

data_set_names = ['mRNA', 'lncRNA', 'sim-mRNA', 'sim-lncRNA']

#Set up plot generator
pg = PlotGenerator()
pg.set_text_options(45, 'right', 0, 'center')

pg.set_text('Mean Base Frequencies', 'Bases', 'Mean Frequency', BASES, None)
pg.bar_plot([mean_pc_base_frequencies, mean_nc_base_frequencies, mean_sim_pc_base_frequencies, mean_sim_nc_base_frequencies], data_set_names)

pg.set_text_options(90, 'right', 0, 'center')

#Bar plots
pg.set_text('Number of Sequences per Sequence Length Range', 'Sequence Length Ranges', 'Number of Sequences', x_axis_labels, None)
pg.bar_plot([pc_bin_sizes, nc_bin_sizes, sim_pc_bin_sizes, sim_nc_bin_sizes], data_set_names)

pg.set_text('Number of Sequences without ORFs', 'Sequence Length Ranges', 'Number of Sequences', x_axis_labels, None)
pg.bar_plot([pc_no_orf_count, nc_no_orf_count, sim_pc_no_orf_count, sim_nc_no_orf_count], data_set_names)

pg.set_text('Number of Sequences of Max ORF Length Equal to or Less than 100', 'Sequence Length Ranges', 'Number of Sequences', x_axis_labels, None)
pg.bar_plot([pc_max_orf_len_less_than_100, nc_max_orf_len_less_than_100, sim_pc_max_orf_len_less_than_100, sim_nc_max_orf_len_less_than_100], data_set_names)

#Box plots
pg.set_axis_options('linear', 10, 'log', 2)

pg.set_text('Length of Longest ORF in RNA Sequences', 'Sequence Length Ranges', 'ORF Length', x_axis_labels, None)
pg.box_plot([pc_max_len_data, nc_max_len_data, sim_pc_max_len_data, sim_nc_max_len_data], data_set_names, True)

pg.set_text('Number of Non-contained ORFs in RNA Sequences', 'Sequence Length Ranges', 'Number of Non-contained ORFs', x_axis_labels, None)
pg.box_plot([pc_max_cnt_data, nc_max_cnt_data, sim_pc_max_cnt_data, sim_nc_max_cnt_data], data_set_names, True)

pg.set_text('Number of Contained ORFs in RNA Sequences', 'Sequence Length Ranges', 'Number of Contained ORFs', x_axis_labels, None)
pg.box_plot([pc_contain_data, nc_contain_data, sim_pc_contain_data, sim_nc_contain_data], data_set_names, True)


# ## Plotting examples
# [boxplot doc](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html)  
# [boxplot demo](https://matplotlib.org/stable/gallery/pyplots/boxplot_demo_pyplot.html)  
# [heatmap examples](https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib) - scroll down!  
