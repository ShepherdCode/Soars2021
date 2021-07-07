import sys
import traceback
import argparse
import re

def assert_imported_GenCode_preprocess():
    return True

'''
Parses an uncompressed GFF3 file from GenCode.
Makes a list of transcript IDs that probably contain an ORF
based on various English-language indicators in the annotation file.
For ease of use within a Python notebook running on CoLab,
this writes a Python file that loads a dict.

Deprecated.
This logic is messy and not 100% accurate.
Instead of this, use GenCodeTools to parse the FASTA file.
The best ORF indicator is a defline with 5'UTR and 3'UTR positions.
'''

class GenCode_Preprocess():
    def __init__(self,debug=False):
        self.debug=debug
        self.prot_incl_tids={}
        self.filename='GenCode_Protein_Include.py'
    def get_filename(self):
        return self.filename
    def get_prot_incl(self):
        return self.prot_incl_tids
    def write_prot_incl_python(self):
        '''Write executable python.'''
        with open (self.filename,'w') as outf:
            print("prot_incl",end="=",file=outf)
            print(self.prot_incl_tids,file=outf)
    def compute_prot_incl(self,gff_fn):
        prot_incl_tids={}
        '''Regular expression.
        Reqires line start with 'chr' to include 'chr22'
        but exclude comment lines etc.
        Requires line not start with 'chrM' to exclude
        mitochondrial genes which use nonstandard codons.
        Requires third field must equal 'transcript'
        in order to exclude 'gene', 'exon' etc.
        '''
        transcript_prefix=re.compile('^chr[^M]+\s.+\stranscript\s')
        EMPTY=''
        with open(gff_fn,'r') as inf:
            for line in inf:
                if not transcript_prefix.match(line):
                    continue
                tid=EMPTY
                rejection=EMPTY
                col9=line.split('\t')[8]
                for pair in col9.split(';'):
                    if pair[:10]=="transcript":
                        if pair[10:14]=="_id=": # transcript_id=
                            taccession=pair[14:]
                            tid=taccession.split('.')[0]
                        elif pair[10:16]=="_type=": #transcript_type
                            if pair[16:]!="protein_coding":
                                rejection="Not protein coding"
                    if pair[:4]=="tag=":
                        if "cds_start_NF" in pair:
                            rejection="No start codon"
                        if "cds_end_NF" in pair:
                            rejection="No stop codon"
                        if "non_ATG_start" in pair:
                            rejection="Noncanonical start codon"
                if rejection != EMPTY:
                    continue  # do not output this one
                if tid == EMPTY:
                    print("DEATH!")
                    break # this should never happen
                prot_incl_tids[tid]=1
                #print(line.rstrip()) # debug only
        self.prot_incl_tids=prot_incl_tids

def args_parse():
    '''GenCode preprocess.'''
    global args
    parser = argparse.ArgumentParser(
        description='GenCode Preprocessor.')
    parser.add_argument(
        'gff',
        help='input annotation (gff)',
        type=str)
    parser.add_argument(
        '--debug',
        help='print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        gff=args.gff
        debug=args.debug
        tool = GenCode_Preprocess(debug)
        sys.stderr.write("Will write to "+tool.get_filename()+"\n")
        tool.compute_prot_incl(gff)
        incl=tool.get_prot_incl()
        tool.write_prot_incl_python()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
