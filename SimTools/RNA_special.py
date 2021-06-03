'''Special Cases'''

def assert_imported_RNA_special():
    return True

class RNA_Special_Cases():
    def __init__(self):
        self.special={}
    def mitochondria(self):
        '''
        Human mitochondrial genes
        can have non-canonical start and stop codons.
        Maybe exclude them from training (but not testing).
        These transcript IDs have field1=chrM
        and gene_type=protein_coding
        in the gencode_annotation.gff file.
        '''
        self.mitochondria=True
        self.special['ENST00000361390']=1
        self.special['ENST00000361453']=1
        self.special['ENST00000361624']=1
        self.special['ENST00000361739']=1
        self.special['ENST00000361851']=1
        self.special['ENST00000361899']=1
        self.special['ENST00000362079']=1
        self.special['ENST00000361227']=1
        self.special['ENST00000361335']=1
        self.special['ENST00000361381']=1
        self.special['ENST00000361567']=1
        self.special['ENST00000361681']=1
        self.special['ENST00000361789']=1
    def is_special(self,TID):
        return TID in self.special
