import numpy as np

'''
Extract K-mer counts or frequencies from RNA sequences.
TO DO: make the Harvester an option on the constructor
so method calls are the same with and without.
'''

class KmerTools():
    def make_kmer_keys(self,K):
        '''
        Make a list of K-mers.
        Returns list of string.
        '''
        shorter_kmers=['']
        for i in range(K):
            longer_kmers=[]
            for mer in shorter_kmers:
                # No support for N or any non-ACGT bases.
                longer_kmers.append(mer+'A')
                longer_kmers.append(mer+'C')
                longer_kmers.append(mer+'G')
                longer_kmers.append(mer+'T')
            shorter_kmers = longer_kmers
        return shorter_kmers
    def make_kmer_dict(self,keys,init=0):
        '''
        Construct a data structure given a list of K-mers.
        The data structure is designed for counting K-mers.
        Returns a dict initialized to (int) zero (default).
        '''
        return dict.fromkeys(keys,init)
    def make_dict_upto_K(self,max_K):
        '''
        Construct a data structure given (int) max_K.
        Returns a dict, values initialized to zeros,
        with one key for every RNA K-mer length 1 to max_K.
        '''
        keys=self.make_kmer_keys(1)
        for k in range(2,max_K+1):
            keys.extend(self.make_kmer_keys(k))
        counts = self.make_kmer_dict(keys)
        return counts
    def update_count_one_K(self,counts,K,rna,tail=False):
        '''
        Given an RNA (string) and K (int),
        count the K-mers in the RNA, and
        update the given data structure called counts.
        Set tail=True only when using Harvester algorithm.
        '''
        L = len(rna)
        padding=" "*(K-1)
        padded=rna+padding
        for i in range(0,L-K+1):
            kmer=padded[i:i+K]
            # This test guards against unexpected bases like N.
            if kmer in counts.keys():
                counts[kmer] += 1
        if tail and K>1:
            # For Harvester algorithm, last letters are special case.
            # Example: for K=3, analyze the last two letters
            # and count two 1-mers and one 2-mer.
            for start_pos in range(L-K+1,L):
                for end_pos in range(start_pos+1,L+1):
                    kmer=rna[start_pos:end_pos]
                    if kmer in counts.keys():
                        counts[kmer] += 1
        return counts
    def update_count_upto_K(self,counts,max_K,sample,tail=False):
        '''
        Repeat update_count_one_K() for every K in 1 to max_K-1.
        '''
        for i in range(1,max_K+1):
            self.update_count_one_K(counts,i,sample,tail)
        return counts
    def harvest_counts_from_K(self,counts,max_K):
        '''
        Implement the Harvester algorithm.
        This is more efficient than update_count_upto_K().
        Given a data structure containing counts for K=max_K only,
        backfill the counts for all smaller values of K.
        Call this after calling update_count_one_K(tail=True)
        on every sequence in your collection with K=max_K.
        '''
        for kmer in counts.keys():
            klen = len(kmer)
            kcnt = counts[kmer]
            if klen==max_K and kcnt>0:
                for i in range(1,klen):
                    prefix = kmer[:i]
                    counts[prefix] += kcnt
        return counts
    def count_to_frequency(self,counts,max_K):
        '''
        Given an data structure of K-mer counts (int),
        return a dict of key:value, type string:double.
        Clients can use list(cnt.values()) to extract just frequencies.
        '''
        freqs = dict.fromkeys(counts.keys(),0.0)
        for k in range(1,max_K+1):
            tot = 0
            for kmer in counts.keys():
                if len(kmer)==k:
                    tot += counts[kmer]
            for kmer in counts.keys():
                if len(kmer)==k:
                    if tot == 0:
                        freqs[kmer] = 0.0
                    else:
                        freqs[kmer] = 1.0*counts[kmer]/tot
        return freqs
    def seqs_to_kmer_freqs(seqs,max_K):
        tool = KmerTools()
        collection = []
        for seq in seqs:
            # Critical to make new dict for each seq,
            # but avoiding it somehow would speed up this loop.
            counts = tool.make_dict_upto_K(max_K)
            # Last param should be True when using Harvester.
            counts = tool.update_count_one_K(counts,max_K,seq,True)
            # Given counts for K=3, Harvester fills in counts for K=1,2.
            counts = tool.harvest_counts_from_K(counts,max_K)
            fdict = tool.count_to_frequency(counts,max_K)
            freqs = list(fdict.values())
            collection.append(freqs)
        return np.asarray(collection)
