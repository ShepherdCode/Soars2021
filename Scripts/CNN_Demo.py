#!/usr/bin/env python
# coding: utf-8

# # CNN Demo
# TO DO: add a CoLab badge
# 
# Can a convolutional neural network (CNN) be trained to distinguish RNA
# from nucleotide composition alone?
# More specifically, can a CNN learn to classify  
# AT-rich sequence with the label "protein coding"
# from GC-rich sequence with the label "non-coding"?
# 
# This demo uses an RNA sequence simulator.
# The simulator strictly follows a frequency histogram with values for A, C, G, T.
# This is a noise-free simulation.
# 
# The CNN is almost as simple as can be.
# It has one trainiable convolution layer (one dimensional) with 8 filters.
# It has one flatten layer simply to reshape the data.
# It has a trainable fully connected (dense) output layer with 1 neuron.
# More sophisticated models would incorporate embedding, pooling, dropout,
# multiple convolution layers, and multiple dense layers.
# 
# The training regime is also simple.
# The model is trained for a fixed number of epochs.
# More sophisticated training would implement early stopping.
# 
# This model minimizes loss at 5 epochs and overfits by 10 epochs.

# ## Computing Environment Setup

# In[19]:


PC_SEQUENCES=2000   # how many protein-coding sequences
NC_SEQUENCES=2000   # how many non-coding sequences
BASES=55            # how long is each sequence
ALPHABET=4          # how many different letters are possible
INPUT_SHAPE_2D = (BASES,ALPHABET,1) # Conv2D needs 3D inputs
INPUT_SHAPE = (BASES,ALPHABET) # Conv1D needs 2D inputs
FILTERS = 8   # how many different patterns the model looks for
WIDTH = 3   # how wide each pattern is, in bases
STRIDE_2D = (1,1)  # For Conv2D how far in each direction
STRIDE = 1 # For Conv1D, how far between pattern matches, in bases
EPOCHS=12  # how many times to train on all the data
SPLITS=4  # SPLITS=3 means train on 2/3 and validate on 1/3 
FOLDS=1  # train the model this many times (must be 1 to SPLITS)


# In[20]:


import sys
try:
    from google.colab import drive
    IN_COLAB = True
    print("On Google CoLab, mount cloud-local file, get our code from GitHub.")
    PATH='/content/drive/'
    #drive.mount(PATH,force_remount=True)
    #drive.mount(PATH)
    DATAPATH=PATH+'My Drive/data/'  # must end in "/"
    import requests
    r = requests.get('https://raw.githubusercontent.com/ShepherdCode/Soars2021/master/SimTools/RNA_gen.py')
    with open('RNA_gen.py', 'w') as f:
        f.write(r.text)  # writes to cloud local, delete the file later?
    from RNA_gen import *
except:
    print("CoLab not working. On my PC, use relative paths.")
    from SimTools.RNA_gen import *
    IN_COLAB = False
    DATAPATH='data/'  # must end in "/"
    sys.path.append("..") # append parent dir in order to use sibling dirs


# In[21]:


from os import listdir
import time # datetime
import csv
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy import stats  # mode

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.utils import shuffle

from keras.models import Sequential
from keras.layers import Dense,Embedding
from keras.layers import Conv1D,Conv2D
from keras.layers import Flatten,MaxPooling1D,MaxPooling2D
from keras.losses import BinaryCrossentropy
# tf.keras.losses.BinaryCrossentropy

import matplotlib.pyplot as plt
from matplotlib import colors
mycmap = colors.ListedColormap(['red','blue'])  # list color for label 0 then 1
np.set_printoptions(precision=2)

if not assert_imported_RNA_gen():
    print("ERROR: Cannot use RNA_gen.")


# ## Data Preparation

# In[22]:


# print(datetime.datetime.now())
t = time.time()
time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(t))


# In[23]:


def make_all_sequences():
    cgen = Collection_Generator()
    sgen = Sequence_Oracle()
    lgen = Length_Oracle()
    lgen.set_mean(BASES)
    cgen.set_seq_oracle(sgen)
    cgen.set_len_oracle(lgen)
    sgen.set_frequencies([4,3,3,4]) # bias AT
    pc_seqs=cgen.get_sequences(PC_SEQUENCES)
    sgen.set_frequencies([3,4,4,3]) # bias GT
    nc_seqs=cgen.get_sequences(NC_SEQUENCES)
    return pc_seqs,nc_seqs
pc_seqs,nc_seqs = make_all_sequences()


# In[24]:


def prepare_for_learning(pcs,ncs):
    samples = nc_seqs + pc_seqs
    NUM_SAMPLES=len(samples)
    X_shape = (NUM_SAMPLES,BASES,ALPHABET)
    Y_shape = (NUM_SAMPLES,1)
    y=np.concatenate((np.zeros(NC_SEQUENCES,dtype=np.int8),
                      np.ones(PC_SEQUENCES,dtype=np.int8)))
    X=np.zeros(X_shape,dtype=np.int8)
    base_to_dim = {'A':0, 'C':1, 'G':2, 'T':3}
    for s in range(0,NUM_SAMPLES):  # TO DO: speed this up by avoiding loops
        sample = samples[s]
        for b in range(0,BASES): # use len(sample) if length varies
            base = sample[b]
            d = base_to_dim[base]   # TO DO: error on non-ACGT
            X[s,b,d]=1
    X,y = shuffle(X,y,random_state=4200)
    return X,y
X,y = prepare_for_learning(pc_seqs,nc_seqs)
print("X shape:",X.shape)
print("y shape:",y.shape)


# ## Model build, train, test

# In[25]:


def make_DNN():
    print("make_DNN")
    print("input shape:",INPUT_SHAPE)
    dnn = Sequential()
    #dnn.add(Embedding(input_dim=4,output_dim=4)) 
    dnn.add(Conv1D(filters=FILTERS,kernel_size=WIDTH,strides=STRIDE,
                  padding="same",input_shape=INPUT_SHAPE))
    # Data shape: [SAMPLES,BASES,FILTERS]
    #dnn.add(MaxPooling1D())
    dnn.add(Flatten())
    # Data shape: [SAMPLES,BASES*FILTERS]
    dnn.add(Dense(1,activation="sigmoid",dtype=np.float32))   
    dnn.compile(optimizer='adam',loss=BinaryCrossentropy(from_logits=False))
    dnn.build(input_shape=INPUT_SHAPE)
    #ln_rate = tf.keras.optimizers.Adam(learning_rate = LN_RATE)
    #bc=tf.keras.losses.BinaryCrossentropy(from_logits=False)
    #model.compile(loss=bc, optimizer=ln_rate, metrics=["accuracy"])
    return dnn
model = make_DNN()
print(model.summary())


# In[26]:


def do_cross_validation(X,y):
    cv_scores = []
    fold=0
    splitter = KFold(n_splits=SPLITS)  # no shuffle
    for train_index,valid_index in splitter.split(X):
        if fold < FOLDS:
            fold += 1
            X_train=X[train_index] # use iloc[] for dataframe
            y_train=y[train_index]
            X_valid=X[valid_index]
            y_valid=y[valid_index]        
            # Avoid continually improving the same model.
            print("MODEL")
            model = model = make_DNN()
            #bestname=MODEL_FILE+".cv."+str(fold)+".best"
            #mycallbacks = [keras.callbacks.ModelCheckpoint(
            #    filepath=bestname, save_best_only=True, 
            #    monitor='val_accuracy', mode='max')]   
            print("FIT")
            start_time=time.time()
            history=model.fit(X_train, y_train, 
                    epochs=EPOCHS, 
                    verbose=1,  # ascii art 
                    # callbacks=mycallbacks,
                    validation_data=(X_valid,y_valid))
            # THE VALIDATION ABOVE IS JUST FOR SHOW
            end_time=time.time()
            elapsed_time=(end_time-start_time)                        
            print("Fold %d, %d epochs, %d sec"%(fold,EPOCHS,elapsed_time))
            pd.DataFrame(history.history).plot(figsize=(8,5))
            plt.grid(True)
            plt.gca().set_ylim(0,1)
            plt.show()
            #best_model=keras.models.load_model(bestname)
            # THE VALIDATION BELOW IS FOR KEEPS
            #scores = best_model.evaluate(X_VALID, Y_VALID, verbose=0)
            #print("%s: %.2f%%" % (best_model.metrics_names[1], scores[1]*100))
            #cv_scores.append(scores[1] * 100)  
    #print()
    #print("%d-way Cross Validation mean %.2f%% (+/- %.2f%%)" % (fold, np.mean(cv_scores), np.std(cv_scores)))


# In[27]:


do_cross_validation(X,y)


# In[27]:




