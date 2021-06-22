import traceback
import argparse
import re
import random

def assert_imported_RNA_describe():
    return True

class ORF_RE():
    '''Regular Expressions'''
    def __init__(self):
        canonical='ATG((?!TAA|TAG|TGA)[ACGT]{3})*(TAA|TAG|TGA)'
        self.canonical_re=re.compile(canonical)
    def get_three_lengths(self,RNA):
        '''Return length of 5'UTR, longest ORF, and 3'UTR.'''
        utr5_len=0
        orf_len=0
        s=self.canonical_re.search(RNA)
        while s is not None:
            this_len=s.end()-s.start()-3 # exclude the stop codon
            if this_len>orf_len:
                orf_len=this_len
                utr5_len=s.start()
            pos = s.start()+1
            s=self.canonical_re.search(RNA,pos)
        if orf_len>0:
            utr3_len=len(RNA)-utr5_len-(orf_len+3)
        else: # no orf? then each UTR is half
            utr5_len=len(RNA)//2
            utr3_len=len(RNA)-utr5_len
        return (utr5_len,orf_len,utr3_len)
    def get_all_orfs(self,RNA):
        orfs=[]
        pos=0
        s=self.canonical_re.search(RNA)
        while s is not None:
            orfs.append(s.group(0))
            pos = s.start()+1
            s=self.canonical_re.search(RNA,pos)
        return orfs

    def get_codon_number(self, RNA):
        '''
        The function get_codon_numbers returns the number of start
        and stop codons present in the parameter.
        It goes through the sequence and takes three letters in frame
        and checks if it is the Start/Stop codon
        '''
        start=stop=0

        for i in range(len(RNA)-2):
            if((RNA[i]=='C') or (RNA[i]=='G')):
                continue

            result = RNA[i]+ RNA[i+1] + RNA[i+2]
            if(result=='ATG'):
                start = start+1
            if((result =='TAA') or (result == 'TAG') or (result == 'TGA')):
                stop = stop+1
        return ([start,stop])


    def get_number_bases(self, RNA):
        '''
        The function below calculates the number of each bases in
        given RNA as parameter. It returns the result as list.
        '''
        a=c=t=g=0
        for i in range(len(RNA)):
            if(RNA[i] == 'A'):
                a+=1
            if(RNA[i] == 'C'):
                c+=1
            if(RNA[i] == 'T'):
                t+=1
            if(RNA[i] == 'G'):
                g+=1

        return([a,c,t,g])

class ORF_counter():
    '''Assume RNA composed of ACGT upper case no N.'''
    def __init__(self,debug=False):
        self.verbose=debug # debugging
        self.set_sequence('')
    def get_max_cds_len(self):
        return self.max_orf_len + 3
    def get_max_orf_len(self):
        return self.max_orf_len
    def count_maximal_orfs(self):
        return self.num_maximal_orfs
    def count_contained_orfs(self):
        return self.num_contained_orfs
    def set_sequence(self,RNA):
        self.max_orf_len=0
        self.num_maximal_orfs=0
        self.num_contained_orfs=0
        self.prev_end=[0,0,0]
        self.prev_start=[0,0,0]
        if self.verbose:
            print("Seq",RNA)
        self.RNA=RNA
        pos=len(RNA)-3
        if pos<3:
            return  # smallest ORF is ATG TAA
        AG = ['A','G']
        while(pos>=0):
            base=RNA[pos]
            plus1=RNA[pos+1]
            plus2=RNA[pos+2]
            if self.verbose:
                print("Pos",pos,"codon",RNA[pos:pos+3])
            if base=='T' and plus1 in AG and plus2 in AG:
                if plus1=='G':
                    if plus2=='A':
                        self.__orf_ends_at__(pos) # TGA
                    if pos>0 and RNA[pos-1]=='A':
                        self.__orf_starts_at__(pos-1) # ATGX
                else: #   TAA, TAG
                    self.__orf_ends_at__(pos)
                pos=pos-3
            elif base=='C':
                pos=pos-3
            else:
                if base=='A' and plus1=='T' and plus2=='G':
                    self.__orf_starts_at__(pos) # ATGX
                pos=pos-1
    def __orf_ends_at__(self,pos):
        # Assume we are working right-to-left (backward)
        # to do: keep track of frame to reduce modulus calls
        frame=pos%3
        if self.verbose:
            print("Frame",frame,"ORF ends at",pos)
        self.prev_start[frame]=0
        self.prev_end[frame]=pos
    def __orf_starts_at__(self,pos):
        # Assume we are working right-to-left (backward)
        # to do: keep track of frame to reduce modulus calls
        frame=pos%3
        if self.verbose:
            print("Frame",frame,"ORF starts at",pos)
        prev_start=self.prev_start[frame]
        prev_end=self.prev_end[frame]
        if prev_end>0:
            if prev_start>0:
                # previously reported ORF wasn't maximal afterall
                self.num_contained_orfs += 1
            else:
                # this orf is maximal (for now)
                self.num_maximal_orfs += 1
            this_len=prev_end-pos
            if this_len>self.max_orf_len:
                self.max_orf_len=this_len
        self.prev_start[frame]=pos

class Random_Base_Oracle():
    def __init__(self,rna_len,debug=False):
        self.rna_len=rna_len
        self.bases=['A','C','G','T']
        self.debug=debug
    def get_one_sequence(self):
        RNA_LEN=self.rna_len
        BASES=self.bases
        between_bases = random.randint(0,RNA_LEN-6)
        utr5_bases = (RNA_LEN - (between_bases + 6)) // 2
        utr3_bases = RNA_LEN - (utr5_bases + (between_bases+6))
        one_seq  = "".join(random.choices(BASES,k=utr5_bases))
        one_seq += 'ATG'
        one_seq += "".join(random.choices(BASES,k=between_bases))
        random_stop = random.choice(['TAA','TAG','TGA']) # random frame
        one_seq += random_stop
        one_seq += "".join(random.choices(BASES,k=utr3_bases))
        if self.debug and len(one_seq) != self.rna_len:
            print("WRONG LENGTH:",len(one_seq),utr5_bases,between_bases,utr3_bases)
        return one_seq
    def get_partitioned_sequences(self,CDS_LEN,goal_per_class):
        pc_seqs=[]
        nc_seqs=[]
        oc = ORF_counter()
        trials = 0
        pc_cnt = 0
        nc_cnt = 0
        while pc_cnt<goal_per_class or nc_cnt<goal_per_class:
            trials += 1
            one_seq=self.get_one_sequence()
            oc.set_sequence(one_seq)
            cds_len = oc.get_max_cds_len()
            if cds_len >= CDS_LEN and pc_cnt<goal_per_class:
                pc_cnt += 1
                pc_seqs.append(one_seq)
            elif cds_len < CDS_LEN and nc_cnt<goal_per_class:
                nc_cnt += 1
                nc_seqs.append(one_seq)
        if self.debug:
            print ("It took %d trials to reach %d per class."%(trials,goal_per_class))
        return pc_seqs,nc_seqs

class RNA_describer():
    '''Assume RNA composed of ACGT upper case no N.'''
    START='ATG'
    STOPS=['TAA','TAG','TGA']
    CODON_LEN=3 # codon length
    def get_orf_length(self,seq):
        '''Count bases within the ORF starting at position 0.
        Return 0 unless first three bases are ATG.
        Count the start codon but not the stop codon (min len 3).
        Only consider in-frame STOP codons.
        Return 0 if no STOP is found.'''
        clen=RNA_describer.CODON_LEN
        if seq[0:clen] == RNA_describer.START:
            if len(seq)>=clen+clen:
                for pos in range(clen,len(seq)-clen+1,clen):
                    if seq[pos:pos+clen] in RNA_describer.STOPS:
                        return pos
        return 0
    def get_longest_orf(self,seq):
        '''Find longest ATG...TAG in any frame.
        Return tuple (offset,len).
        Offset starts at zero.
        Length includes ATG and TAG (min 6).
        Returns (0,0) if none found.'''
        clen=RNA_describer.CODON_LEN
        rlen = len(seq)
        longest_orf=(0,0) # offset, len
        for one_pos in range(0,rlen):
            # TO DO: an optimization for seq like START,START,STOP.
            # No need to count length at second START.
            if seq[one_pos:one_pos+clen] == RNA_describer.START:
                one_len = self.get_orf_length(seq[one_pos:])
                if one_len > longest_orf[1]:
                    longest_orf=(one_pos,one_len)
        return longest_orf
    def get_orf_lengths(self,seqs):
        '''Given list of zero to many sequences.
        Return list of largest ORF length per sequence.'''
        lens = []
        for one_seq in seqs:
            # orf=(offset,length)
            one_orf = self.get_longest_orf(one_seq)
            one_len = one_orf[1]
            lens.append(one_len)
        return lens
    def get_three_lengths(self,seqs):
        '''Given list of zero to many sequences.
        Find and use the longest ORF per sequence.
        Return list of three lengths per sequence, like this:
        [ (5'UTR, ORF, 3'UTR) ].
        This supports statistical characterization of a data set.
        For sequences lacking an ORF, return (len/2,0,len/2).'''
        bounds = []
        for one_seq in seqs:
            # orf=(offset,length)
            one_orf = self.get_longest_orf(one_seq)
            seq_length=len(one_seq)
            orf_length=one_orf[1]
            if orf_length==0:
                one_bound = ((seq_length+1)//2,0,seq_length//2)
            else:
                # 5'UTR ends where ORF begins
                utr5_length=one_orf[0]
                # 3'UTR begins after the stop codon
                # (which was not counted as part of the ORF)
                utr3_length=seq_length-utr5_length-orf_length-3
                one_bound = (utr5_length,orf_length,utr3_length)
            bounds.append(one_bound)
        return bounds

class FASTA_tool():
    def __init__(self,infile,debug=False):
        self.fn=infile
        self.verbose=debug
    def show_all(self):
        prev=None
        seq=None
        oc=ORF_counter()
        with open(self.fn,'r') as inf:
            for line in inf:
                sline=line.rstrip()
                if sline[0]=='>':
                    if seq is not None:
                        rna_len=len(seq)
                        oc.set_sequence(seq)
                        max_orf_len=oc.get_max_orf_len()
                        print(prev,rna_len,max_orf_len)
                    prev=sline[1:].split('|')[0]
                    seq=""
                else:
                    seq += sline
            if seq is not None:
                rna_len=len(seq)
                oc.set_sequence(seq)
                max_orf_len=oc.get_max_orf_len()
                print(prev,rna_len,max_orf_len)


def args_parse():
    '''RNA describer.'''
    global args
    parser = argparse.ArgumentParser(
        description='RNA Describer.')
    parser.add_argument(
        'infile',
        help='input filename (fasta)',
        type=str)
    parser.add_argument(
        '--debug',
        help='print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    try:
        args_parse()
        infile=args.infile
        debug=args.debug
        tool = FASTA_tool(infile,debug)
        tool.show_all()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
