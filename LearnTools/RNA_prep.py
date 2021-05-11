from sklearn.utils import shuffle
import numpy as np

def prepare_inputs_len_x_alphabet(
pc_seqs,nc_seqs,alphabet_size,with_shuffle=True):
    samples = nc_seqs + pc_seqs
    seq_len=len(nc_seqs[0]) # TO DO: error on empty or non-uniform
    num_samples=len(samples)
    X_shape = (num_samples,seq_len,alphabet_size)
    Y_shape = (num_samples,1) # one label per sequence
    y=np.concatenate((np.zeros(len(nc_seqs),dtype=np.int8),
                      np.ones(len(pc_seqs),dtype=np.int8)))
    X=np.zeros(X_shape,dtype=np.int8)
    base_to_dim = {'A':0, 'C':1, 'G':2, 'T':3}
    for s in range(0,num_samples):  # TO DO: speed this up by avoiding loops
        sample = samples[s]
        for b in range(0,seq_len):
            base = sample[b]
            d = base_to_dim[base]   # TO DO: error on non-ACGT
            X[s,b,d]=1
    if with_shuffle:
        X,y = shuffle(X,y,random_state=4200)
    return X,y
