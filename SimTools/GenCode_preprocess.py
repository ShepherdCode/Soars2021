import traceback
import argparse
import re

def assert_imported_GenCode_preprocess():
    return True

class GenCode_Preprocess():
    def __init__(self,debug=False):
        self.debug=debug
    def get_prot_incl(self,gff_fn):
        prot_incl_tids=[]
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
                if rejection != EMPTY:
                    continue  # do not output this one
                if tid == EMPTY:
                    print("DEATH!")
                    break # this should never happen
                prot_incl_tids.append(tid)
                print(line.rstrip()) # debug only
        return prot_incl_tids

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
        incl=tool.get_prot_incl(gff)
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
