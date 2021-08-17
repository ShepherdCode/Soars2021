import numpy as np
from sklearn.utils import shuffle

class DataPrep():
    def combine_pos_and_neg(self,seqs1,seqs0):
        len1=len(seqs1)
        len0=len(seqs0)
        if False:  # high ram footprint high due to variable length seqs
            L1=np.ones(len1,dtype=np.int8)
            L0=np.zeros(len0,dtype=np.int8)
            S1 = np.asarray(seqs1)
            S0 = np.asarray(seqs0)
            all_labels = np.concatenate((L1,L0))
            all_seqs = np.concatenate((S1,S0))
        else:  # numpy free saves memory
            L1=[1]*len1
            L0=[0]*len0
            S1=list(seqs1)
            S0=list(seqs0)
            all_seqs=S1+S0
            all_labels=L1+L0
        X,y = shuffle(all_seqs,all_labels) # sklearn.utils.shuffle
        return X,y
