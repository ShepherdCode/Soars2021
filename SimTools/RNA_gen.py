import traceback
import argparse
# By default, the random library uses one
# singleton instance shared by all client code.
import random

def assert_imported_RNA_gen():
    return True

class Length_Oracle():
    '''Generate one sequence length.
    Always returns the same number.
    This is intended as a base class.'''
    def __init__(self):
        self.set_mean(50)
    def set_mean(self,mean):
        self.mean=mean
    def get_length(self):
        return self.mean

class Sequence_Oracle():
    '''Generate one RNA sequence.'''
    def __init__(self):
        self.set_sequences()
        self.set_frequencies()
        self.check_params()
    def set_sequences(self,seqs=['A','C','G','T']):
        self.seqs=seqs
    def set_frequencies(self,freqs=[1,1,1,1]):
        self.freqs=freqs
    def check_params(self):
        seqs=self.seqs
        freqs=self.freqs
        if len(seqs)==0 or len(freqs)==0:
            print("WARN: missing seqs or freqs")
            return False
        if not len(seqs) == len(freqs):
            print("WARN: len(seqs)!=len(freqs)")
            return False
        # TO DO: check for valid freqs
    def get_sequence(self,len):
        rnd = random.choices(self.seqs,
            weights=self.freqs,k=len)
        seq=''.join(rnd)
        return seq

class Transcript_Oracle(Sequence_Oracle):
    def __init__(self):
        self.orf_oracle = Sequence_Oracle()
        # TO DO: change next line to use codons
        self.orf_oracle.set_frequencies([0,0,0,1])
        self.utr_portion=3  # one third
        super().__init__()
    def get_sequence(self,len):
        utr_len = len//self.utr_portion
        orf_len = len - utr_len - utr_len
        utr5 = super().get_sequence(utr_len)
        utr3 = super().get_sequence(utr_len)
        orf = self.orf_oracle.get_sequence(orf_len)
        return utr5 + orf + utr3

class Collection_Generator():
    '''Generate one file of simulated RNA.'''
    def __init__(self,debug=False):
        self.debug=debug
        self.filename="Generated_RNA.fasta"
        self.set_seq_oracle(Sequence_Oracle())
        self.set_len_oracle(Length_Oracle())
    def set_reproducible(self,same):
        if same:  # reproducible behavior
            REPRODUCIBLE_SEED=42
            random.seed(REPRODUCIBLE_SEED)
        else:  # default pseudo-random behavior
            random.seed(None)
    def set_filename(self,fn):
        self.filename=fn
    def set_seq_oracle(self,so):
        self.sequence_oracle = so
    def set_len_oracle(self,lo):
        self.length_oracle=lo
    def get_seq_oracle(self):
        return self.sequence_oracle
    def get_len_oracle(self):
        return self.length_oracle
    def get_sequences(self,seqs=1):
        lo = self.length_oracle
        so = self.sequence_oracle
        collection=[]
        for num in range(0,seqs):
            len = lo.get_length()
            seq = so.get_sequence(len)
            collection.append(seq)
        return collection
    def write_fasta(self,seqs=10):
        fn = self.filename
        lo = self.length_oracle
        so = self.sequence_oracle
        with open(fn, 'w') as outfa:
            for num in range(0,seqs):
                id="random_%07d"%num
                outfa.write(">%s\n"%id)
                len = lo.get_length()
                seq = so.get_sequence(len)
                outfa.write(seq+"\n")

def args_parse():
    '''RNA simulator.'''
    global args
    parser = argparse.ArgumentParser(
        description='RNA Simulator.')
    parser.add_argument(
        'numlines',
        help='output file size (10)',
        type=int)
    parser.add_argument(
        'outfile',
        help='output filename (fasta)',
        type=str)
    parser.add_argument(
        '--transcript',
        help='create UTR/ORF/UTR.',
        action='store_true')
    parser.add_argument(
        '--debug',
        help='print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        numlines=args.numlines
        outfile=args.outfile
        debug=args.debug
        gen = Collection_Generator(debug)
        if args.transcript:
            gen.set_seq_oracle(Transcript_Oracle())
        gen.set_filename(outfile)
        gen.write_fasta(numlines)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
