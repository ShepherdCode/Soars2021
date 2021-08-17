import pytest
from ORF_eliminator import ORF_eliminator
from RNA_describe import RNA_describer
from RNA_describe import ORF_counter
from RNA_describe import ORF_RE

class Test_ORF_eliminator:
    '''
    Test to see if ORF_eliminator recreates the sequence
    without ORF!
    In rna, clearly we have a ORF of length 6.
    We will check if goes down to 0 using get_coordinate function.
    
    Also, we will check the length to see if the length of new
    RNA is same or not. 
'''
    def test_method1(self):
        rna = 'ATGAAATAG'
        length_rna = len(rna)
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF()
        length_new_rna = len(newRNA)
        coordinates =data.get_coordinates()
        length_of_orfs = coordinates[1]
        ans = 0
        msg ='Removed!'
        msg2 = 'Both RNA are of same length!'
        assert length_of_orfs == ans , msg
        assert length_rna == length_new_rna ,msg2 
        '''
        This time, lets make the RNA with just codon TAG
        rna = 'TAGTAGTAGTAGTAG'

'''
    def test_method2(self):
        rna = 'TAGTAGTAGTAGTAG'
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF()
        coordinate =data.get_coordinates()
        length = coordinate[1]
        ans = 0
        msg ='Removed!'
        assert length == ans , msg
        
        length_rna = len(rna)
        length_new_rna = len(newRNA)
        msg2 = 'Both RNA are of same length!'
        assert length_rna == length_new_rna ,msg2 
    '''
        Lets use letter T all the time to see if the
        program fails. 
    '''
    
    def test_method3(self):
        rna = 'TTTTTTTTTTTTTTT'
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF()
        coordinate =data.get_coordinates()
        length = coordinate[1]
        ans = 0
        msg ='Removed!'
        assert length == ans , msg
        length_rna = len(rna)
        length_new_rna = len(newRNA)
        msg2 = 'Both RNA are of same length!'
        assert length_rna == length_new_rna ,msg2 
        
    def test_method4(self):
        rna = 'ATGAAATAG'
        length_rna = len(rna)
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF2()
        length_new_rna = len(newRNA)
        coordinates =data.get_coordinates()
        length_of_orfs = coordinates[1]
        ans = 0
        msg ='Removed!'
        msg2 = 'Both RNA are of same length!'
        assert length_of_orfs == ans , msg
        assert length_rna == length_new_rna ,msg2 
        '''
        This time, lets make the RNA with just codon TAG
        rna = 'TAGTAGTAGTAGTAG'

'''
    def test_method5(self):
        rna = 'TAGTAGTAGTAGTAG'
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF2()
        coordinate =data.get_coordinates()
        length = coordinate[1]
        ans = 0
        msg ='Removed!'
        assert length == ans , msg
        
        length_rna = len(rna)
        length_new_rna = len(newRNA)
        msg2 = 'Both RNA are of same length!'
        assert length_rna == length_new_rna ,msg2 
    '''
        Lets use letter T all the time to see if the
        program fails. 
    '''
    
    def test_method6(self):
        rna = 'TTTTTTTTTTTTTTT'
        data = ORF_eliminator(rna)
        newRNA = data.eliminate_ORF2()
        coordinate =data.get_coordinates()
        length = coordinate[1]
        ans = 0
        msg ='Removed!'
        assert length == ans , msg
        length_rna = len(rna)
        length_new_rna = len(newRNA)
        msg2 = 'Both RNA are of same length!'
        assert length_rna == length_new_rna ,msg2 
        