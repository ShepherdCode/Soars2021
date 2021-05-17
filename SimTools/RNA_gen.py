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
    Intended as a base class.'''
    def __init__(self):
        self.mean=50 # arbitrary default
        self.stddev=0
    def set_mean(self,mean):
        self.mean=mean
    def set_stddev(self,std):
        self.stddev=std
    def get_length(self):
        value=self.mean
        if self.stddev>0:
            value=random.gauss(self.mean,self.stddev)
            value=int(value+0.5)
        return value

class Sequence_Oracle():
    '''Generate one RNA sequence.
    Use weighted random selection of given K-mers.
    Intended as a base class.'''
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
    def get_sequence(self,length):
        rnd = random.choices(self.seqs,
            weights=self.freqs,k=length)
        seq=''.join(rnd)
        return seq

class Transcript_Oracle(Sequence_Oracle):
    '''Return a sequence containing an ORF.
    The ORF is centered along the transcript.
    The start is always ATG but the stop is randomly chosen.
    The codons are uniform random but without stops.'''
    def __init__(self,debug=False):
        super().__init__()
        self.debug = debug
        self.orf_len=None # by default, use portion of transcript
        self.orf_portion=2  # ORF len = transcript len / 2
        self.utr5_portion=2 # UTR5 len = total UTR / 2
        self.min_orf_len=9  # START+CODON+STOP
        self.var=1/6  # orf len stddev/mean
        self.codon_size=3
        self.START = 'ATG'
        self.STOPS = ['TAA','TAG','TGA']
        codons = Transcript_Oracle.make_codons('ACGT')
        codons = Transcript_Oracle.remove_stops(codons,self.STOPS)
        freqs = [1]*len(codons)
        # subclasses could use change these
        self.orflen_oracle = Length_Oracle()
        self.codons_oracle = Sequence_Oracle()
        self.codons_oracle.set_sequences(codons)
        self.codons_oracle.set_frequencies(freqs)
    def make_codons(bases):
        codons=[]
        num_bases=len(bases)
        for i in range(0,num_bases):
            for j in range(0,num_bases):
                for k in range(0,num_bases):
                    codon=bases[i]+bases[j]+bases[k]
                    codons.append(codon)
        return codons
    def remove_stops(codons,stops):
        for stop in stops:
            codons.remove(stop)
        return codons
    def set_orf_len_mean(self,value):
        self.orf_len=(value)
    def get_sequence(self,length):
        '''Generates 5'UTR + ORF + 3'UTR.
        Both UTR contain random sequence.
        ORF structure is START+CODON(S)+STOP.'''
        orf = self.__get_orf(length)
        utr5,utr3 = self.__get_utr(length,len(orf))
        if self.debug:
            orf = orf.lower()  # for visual debugging
        return utr5 + orf + utr3
    def __get_orf(self,tlen):
        if self.orf_len is None:
            orflen_target = tlen//self.orf_portion
        else:
            orflen_target = self.orf_len # user settable
        self.orflen_oracle.set_mean(orflen_target)
        self.orflen_oracle.set_stddev(orflen_target*self.var)
        orflen_actual = self.orflen_oracle.get_length()
        orflen_actual = min(tlen,orflen_actual)
        orflen_actual -= orflen_actual%self.codon_size
        if orflen_actual<self.min_orf_len:
            # Too short. Just return random sequence.
            return super.get_sequence(tlen)
        num_codons = orflen_actual//self.codon_size
        num_codons -= 2 # START and STOP are automatic
        orf = self.START
        orf += self.codons_oracle.get_sequence(num_codons)
        orf += random.choice(self.STOPS)
        return orf
    def __get_utr(self,tlen,orflen):
        utr5 = ''
        utr3 = ''
        utr5_len = (tlen-orflen)//self.utr5_portion
        utr3_len = tlen - orflen - utr5_len
        if utr5_len > 0:
            utr5 = super().get_sequence(utr5_len)
        if utr3_len > 0:
            utr3 = super().get_sequence(utr3_len)
        return utr5,utr3

class Collection_Generator():
    '''Generate one collection of simulated RNA.
    Optionally write the collection to a FASTA file.'''
    def __init__(self,debug=False):
        self.debug=debug  # needed?
        self.set_seq_oracle(Sequence_Oracle())
        self.set_len_oracle(Length_Oracle())
    def set_reproducible(self,same):
        if same:  # reproducible behavior
            REPRODUCIBLE_SEED=42
            random.seed(REPRODUCIBLE_SEED)
        else:  # default pseudo-random behavior
            random.seed(None)
    def set_seq_oracle(self,so):
        self.sequence_oracle = so
    def set_len_oracle(self,lo):
        self.length_oracle=lo
    def get_seq_oracle(self):
        return self.sequence_oracle
    def get_len_oracle(self):
        return self.length_oracle
    def get_sequences(self,seqs=1):
        '''Return a list of simulated RNA strings.'''
        lo = self.length_oracle
        so = self.sequence_oracle
        collection=[]
        for num in range(0,seqs):
            len = lo.get_length()
            seq = so.get_sequence(len)
            collection.append(seq)
        return collection
    def write_fasta(self,filename,seqs=10):
        '''Write a file of simulated RNA strings.'''
        lo = self.length_oracle
        so = self.sequence_oracle
        with open(filename, 'w') as outfa:
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
        help='desired sequence count',
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
            gen.set_seq_oracle(Transcript_Oracle(debug))
        gen.write_fasta(outfile,numlines)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
