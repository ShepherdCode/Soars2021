import re  # regular expressions
import gzip
import pandas as pd

class GenCodeLoader():
    def __init__(self):
        self.pattern5=re.compile('.*UTR5:')
        self.pattern3=re.compile('.*UTR3:')
        self.check_list = None
        self.check_utr = False
    def set_label(self,label):
        self.label=label
    def set_check_list(self,check_list):
        self.check_list=check_list
    def set_check_utr(self,check_utr):
        self.check_utr=check_utr
    def __save_previous(self,one_def,one_seq):
        if one_def is None:
            return
        if self.check_utr:
            if self.pattern5.match(one_def) is None:
                return
            if self.pattern3.match(one_def) is None:
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
