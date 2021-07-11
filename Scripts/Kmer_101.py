#!/usr/bin/env python
# coding: utf-8

# # K-mer 
# Move code into SimTools library.

# In[1]:


import time
def show_time():
    t = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t)))
show_time()


# In[2]:


PC_SEQUENCES=32000
NC_SEQUENCES=32000
RNA_LEN=32
CDS_LEN=16


# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[4]:


import sys
IN_COLAB = False
try:
    from google.colab import drive
    IN_COLAB = True
except:
    pass
if IN_COLAB:
    print("On Google CoLab, mount cloud-local file, get our code from GitHub.")
    PATH='/content/drive/'
    #drive.mount(PATH,force_remount=True)  # hardly ever need this
    #drive.mount(PATH)    # Google will require login credentials
    DATAPATH=PATH+'My Drive/data/'  # must end in "/"
    import requests
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_describe.py')
    with open('RNA_describe.py', 'w') as f:
        f.write(r.text)  
    from RNA_describe import Random_Base_Oracle
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/KmerTools.py')
    with open('KmerTools.py', 'w') as f:
        f.write(r.text)  
    from KmerTools import KmerTools
else:
        print("CoLab not working. On my PC, use relative paths.")
        DATAPATH='data/'  # must end in "/"
        sys.path.append("..") # append parent dir in order to use sibling dirs
        from SimTools.RNA_describe import Random_Base_Oracle
        from SimTools.KmerTools import KmerTools
MODELPATH="BestModel"  # saved on cloud instance and lost after logout
#MODELPATH=DATAPATH+MODELPATH  # saved on Google Drive but requires login


# ## Demo

# ### Use Harvester algorithm

# In[5]:


rbo=Random_Base_Oracle(RNA_LEN,True)
pc_all,nc_all = rbo.get_partitioned_sequences(CDS_LEN,10) # just testing
pc_all,nc_all = rbo.get_partitioned_sequences(CDS_LEN,PC_SEQUENCES)
print("Use",len(pc_all),"PC seqs")
print("Use",len(nc_all),"NC seqs")


# In[6]:


MAX_K = 3
tool = KmerTools()
pc_counts = tool.make_dict_upto_K(MAX_K)
for sample in pc_all:
    tool.update_count_one_K(pc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(pc_counts,MAX_K)
print("PC counts:\n",pc_counts)
pc_freqs = tool.count_to_frequency(pc_counts,MAX_K)
print ("Frequency:\n",pc_freqs)


# In[7]:


nc_counts = tool.make_dict_upto_K(MAX_K)
for sample in nc_all:
    tool.update_count_one_K(nc_counts,MAX_K,sample,True)
tool.harvest_counts_from_K(nc_counts,MAX_K)
print("NC counts:\n",nc_counts)
nc_freqs = tool.count_to_frequency(nc_counts,MAX_K)
print ("Frequency:\n",nc_freqs)


# In[ ]:




