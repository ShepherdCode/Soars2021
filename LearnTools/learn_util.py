from sklearn.utils import shuffle
def prepare_for_learning(pcs,ncs,with_shuffle=True):
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
    if with_shuffle:
        X,y = shuffle(X,y,random_state=4200)
    return X,y
