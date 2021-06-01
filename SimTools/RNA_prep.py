from sklearn.utils import shuffle
import numpy as np

def assert_imported_RNA_prep():
    return True

def prepare_inputs_len_x_alphabet(pc_seqs,nc_seqs,alphabet_size,with_shuffle=True):
    # TO DO: eliminate alphabet_size and hard code the 4?
    samples = nc_seqs + pc_seqs
    num_samples=len(samples)
    # Assume all sequences have same length
    seq_len=len(samples[0]) # TO DO: allow non-uniform size
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

def get_codon_index():
    codon_to_index={}
    index = 0
    for b1 in 'ACGT':
        for b2 in 'ACGT':
            for b3 in 'ACGT':
                codon=b1+b2+b3
                codon_to_index[codon]=index
                index += 1
    return codon_to_index

def prepare_inputs_codon_frequency(pc_seqs,nc_seqs,with_shuffle=True):
    alphabet_size=4
    codon_size=3
    codon_space = alphabet_size**codon_size # 4^3=64
    samples = nc_seqs + pc_seqs
    num_samples=len(samples)
    X_shape = (num_samples,codon_space)
    Y_shape = (num_samples,1) # one label per sequence
    y=np.concatenate((np.zeros(len(nc_seqs),dtype=np.int8),
                      np.ones(len(pc_seqs),dtype=np.int8)))
    X=np.zeros(X_shape,dtype=np.float32)
    codon_to_index=get_codon_index()
    for s in range(0,num_samples):
        sample = samples[s]
        thousandth = 1/len(sample)
        seq_len=len(sample)
        for b in range(0,seq_len-codon_size+1):
            codon = sample[b:b+codon_size]
            index = codon_to_index[codon]
            X[s,index]+=thousandth
    if with_shuffle:
        X,y = shuffle(X,y,random_state=4200)
    return X,y
