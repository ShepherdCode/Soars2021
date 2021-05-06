import os, sys, traceback, argparse

class Length_Oracle():
    '''Generate one sequence length.'''
    def get_length(self):
        return 50

class Sequence_Oracle():
    '''Generate one RNA sequence.'''
    def __init__(self,debug=False):
        self.debug=debug
    def get_sequence(self,len):
        seq = "A"*len
        return seq

class File_Generator():
    '''Generate one file of simulated RNA.'''
    def __init__(self,debug=False):
        self.debug=debug
        self.filename="Generated_RNA.fasta"
        self.seq_oracle=Sequence_Oracle()
        self.len_oracle=Length_Oracle()
    def set_filename(self,fn):
        self.filename=fn
    def write_file(self,seqs=10):
        fn = self.filename
        with open(fn, 'w') as outfa:
            for num in range(0,seqs):
                outfa.write(">sequence_"+str(num)+"\n")
                len = self.len_oracle.get_length()
                seq = self.seq_oracle.get_sequence(len)
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
        '--debug',
        help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        numlines=args.numlines
        outfile=args.outfile
        debug=args.debug
        gen = File_Generator(debug)
        gen.set_filename(outfile)
        gen.write_file(numlines)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
