import traceback
import argparse
import re  # regular expressions
import gzip
import pandas as pd

'''
Load RNA sequence into memory.
Reads a FASTA.gz file from GeneCode.
Parses the transcript id (TID) from the FASTA defline.
Returns a Pandas dataframe with columnts tid, class, sequence, seqlen.
Typical input files from (https://www.gencodegenes.org/)
- gencode.v38.lncRNA_transcripts.fa.gz
- gencode.v38.pc_transcripts.fa.gz
'''

class GenCodeLoader():
    def __init__(self):
        self.pattern5=re.compile('.*UTR5:')
        self.pattern3=re.compile('.*UTR3:')
        self.check_list = None
        self.check_utr = False
        self.min_size = None
        self.max_size = None
    def set_label(self,label):
        '''
        Set one label used for subsequent sequences.
        The value gets stored in the 'class' field.
        Usually use 1 for protein-coding and 0 for non-coding.
        '''
        self.label=label
    def set_check_list(self,check_list):
        '''
        Optionally provide a TID include list. Others are excluded.
        The parameter, type list, is used with pythin 'in' operator.
        '''
        self.check_list=check_list
    def set_check_utr(self,check_utr):
        '''
        Optionally require UTR. Equivalent to requiring an ORF.
        Include only deflines that specify 5'UTR and 3'UTR positions.
        (GenCode does have mRNA transcripts that lack an ORF!)
        Set this to false when loading non-coding RNA.
        The parameter is type boolean.
        '''
        self.check_utr=check_utr
    def set_check_size(self,min,max):
        self.min_size = min
        self.max_size = max
    def __save_previous(self,one_def,one_seq):
        '''
        For internal use only.
        FASTA sequence records are multi-line starting with a defline.
        This is called just before parsing a new defline
        to optionally save the previously parsed sequence record.
        '''
        if one_def is None:
            return
        if self.check_utr:
            if self.pattern5.match(one_def) is None:
                return
            if self.pattern3.match(one_def) is None:
                return
        seq_len = len(one_seq)
        if self.min_size is not None and seq_len < self.min_size:
            return
        if self.max_size is not None and seq_len > self.max_size:
            return
        VERSION = '.'
        one_id = one_def[1:].split(VERSION)[0]
        if self.check_list is not None:
            if one_id not in self.check_list:
                return
        self.labels.append(self.label)
        self.seqs.append(one_seq)
        self.lens.append(len(one_seq))
        self.ids.append(one_id)
    def load_file(self,filename):
        '''
        Parse the given file and return a data structure.
        Given file assumed GenCode FASTA file.
        Returns a Pandas dataframe with four fields.
        '''
        self.labels=[]  # usually 1 for protein-coding or 0 for non-coding
        self.seqs=[]    # usually strings of ACGT
        self.lens=[]    # sequence length
        self.ids=[]     # GenCode transcript ID, always starts ENST, excludes version
        DEFLINE='>'  # start of line with ids in a FASTA FILE
        EMPTY=''
        one_def = None
        one_seq = ''
        with gzip.open (filename,'rt') as infile:
            for line in infile:
                if line[0]==DEFLINE:
                    self.__save_previous(one_def,one_seq)
                    one_def=line
                    one_seq = EMPTY
                else:
                    # Continue loading sequence lines till next defline.
                    additional = line.rstrip()
                    one_seq = one_seq + additional
            # Don't forget to save the last sequence after end-of-file.
            self.__save_previous(one_def,one_seq)
        df1=pd.DataFrame(self.ids,columns=['tid'])
        df2=pd.DataFrame(self.labels,columns=['class'])
        df3=pd.DataFrame(self.seqs,columns=['sequence'])
        df4=pd.DataFrame(self.lens,columns=['seqlen'])
        df=pd.concat((df1,df2,df3,df4),axis=1)
        return df

def args_parse():
    '''GenCode FASTA parser.'''
    global args
    parser = argparse.ArgumentParser(
        description='GenCode FASTA parser.')
    parser.add_argument(
        'infile',
        help='input filename (gzip fasta)',
        type=str)
    parser.add_argument(
        '--orf',
        help='require UTR/ORF/UTR in defline',
        action='store_true')
    parser.add_argument(
        '--debug',
        help='print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        parser = GenCodeLoader()
        parser.set_label(1)
        if args.orf:
            parser.set_check_utr(True)
        df = parser.load_file(args.infile)
        print(df)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
