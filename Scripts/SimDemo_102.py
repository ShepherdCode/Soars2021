#!/usr/bin/env python
# coding: utf-8

# # RNA Simulator from NASA Spring 2021
# This simulator was developed by Keagen Thomson during his NASA fellowship with John Little. 
# 
# The original ORF generator was at:[ShepherdCode/SheperdML/Nasa2021/GenFASTA.py](https://github.com/ShepherdCode/ShepherdML/blob/master/Nasa2021/GenFASTA.py)
# 
# The companion ORF finder was at:[ShepherdCode/SheperdML/Nasa2021/PyORF.py]()
# 
# TO DO: Transform this code from a notebook into a SimTools module.

# In[1]:


import random
import itertools
#import os
import sys

sys.path.append("..")
from SimTools.RNA_describe import ORF_counter


# In[2]:


Pairs = ("A","C", "G", "T")
Codons = itertools.product(Pairs, repeat=3)
all_codons = ["".join(codon) for codon in Codons]
start = "ATG"
stop = ("TAA", "TAG", "TGA") 
no_stops = [codon for codon in all_codons if codon not in stop]


# Prepends a random A, G, C, or T for the input frame e.g. prepends 1 character if in frame 2, and 2 if in frame 3
# 

# In[3]:


def shift_frame(input_seq,frame = 2):
  output = input_seq
  if frame in (1,2,3):
    for i in range(frame-1):
      output.insert(0, random.choice(("A","G","C","T")))
    return output
  else:
    raise ValueError("Frame Must Be 1, 2 or 3. Frame Entered: " +frame)


# Places random codons length times, and uses choices from no_stops if coding, and all_codons if not coding

# In[4]:


def codon_placer(length, coding = True):
  lno_stops = no_stops
  lall_codons = all_codons

  if coding == True:
    return random.choices(lno_stops, k=length)
  else:
    return random.choices(lall_codons,k=length)


# Returns random indexes a start and stop codon should be placed arbitrarily chooses a variation within 1/3rd the 2nd and 4th Quintile

# In[5]:


def get_index_placement(total_codons):
  quintile = total_codons // 5
  variation = quintile // 3
  start_placement = quintile + random.randint(-variation, variation)
  stop_placement = (quintile*4) + random.randint(-variation, variation)
  return start_placement, stop_placement


# Generates a random (hypothesized) coding or non-coding sequence of length characters in frame  
# length: Number of characters the sequence should be  
# coding: Whether or not the sequence should have stop codons placed, or randomly generated (True, False respectively)  
# frame: The frame the sequence should be in which determines how many bases are appended to the sequences start. Must be 1, 2 or 3

# In[6]:


def generate_seq(length, coding = False, frame = 1):
  codons_to_place = (length//3) + 1

  if coding and frame in (1,2,3):
    start_index, stop_index = get_index_placement(codons_to_place)
    UTR_5_len = start_index-1
    orf_length = stop_index-start_index - 2
    UTR_3_len = codons_to_place - stop_index + 1
      
    UTR_5 = codon_placer(UTR_5_len, False)
    sequence_orf = codon_placer(orf_length, True)
    sequence_orf.insert(0, start)
    sequence_orf.append(random.choice(stop))
    UTR_3 = codon_placer(UTR_3_len, False)

    UTR_5.extend(sequence_orf)
    UTR_5.extend(UTR_3)
    output = shift_frame(UTR_5, frame)
    output = ''.join(output)
    # Unclear what last return pair represented.
    #return output[0:length], coding, frame, (start_index, stop_index)
    return output[0:length], coding, frame

  elif not coding and frame in (1,2,3):
    placed_codons = codon_placer(codons_to_place, coding)
    output = shift_frame(placed_codons, frame)
    output = ''.join(output)
    return output[0:length] , coding, frame
  else:
    raise ValueError("Frame must be 1, 2 or 3")


# Try using the code above.

# In[7]:


PC=True
NC=False
RNA_LEN=15


# In[8]:


seq1=generate_seq(RNA_LEN,PC,2)
print(seq1)


# In[9]:


RNA_LEN=1000
PC_SAMPLES=1000
NC_SAMPLES=1000
pc_all=[]
nc_all=[]
frame=1
for i in range(0,PC_SAMPLES):
    rec=generate_seq(RNA_LEN,NC)
    seq=rec[0]
    nc_all.append(seq)
    rec=generate_seq(RNA_LEN,PC,frame)
    seq=rec[0]
    pc_all.append(seq)
    frame += 1
    if frame>3:
        frame=1
print("PC count",len(pc_all))
print(pc_all[0])
print("NC count",len(nc_all))
print(nc_all[0])    


# Evaluate the results above

# In[10]:


oc = ORF_counter()
print("Simulated sequences:")
print("PC seqs")
oc.describe_sequences(pc_all)
print("NC seqs")
oc.describe_sequences(nc_all)

