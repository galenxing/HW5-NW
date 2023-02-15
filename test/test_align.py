# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    nw.align(seq1,seq2)

    assert nw._align_matrix.shape == (len(seq1)+1 , len(seq2)+1)
    assert nw._gapA_matrix.shape == (len(seq1)+1 , len(seq2)+1)
    assert nw._gapB_matrix.shape == (len(seq1)+1 , len(seq2)+1)
    
    assert nw._gapA_matrix[2:,1].all() == True
    assert nw._gapB_matrix[1,2:].all() == True
    
    assert nw._align_matrix[len(seq1),len(seq2)] == 4
    assert nw._align_matrix[0,0] == 0
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    score, s1, s2 = nw.align(seq3,seq4)
  
    assert score == 17
    assert s1 == 'MAVHQLIRRP'
    assert s2 == 'M---QLIRHP'
    


    
    




