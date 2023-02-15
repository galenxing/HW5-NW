# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)

    gg,_,_ = nw.align(hs_seq,gg_seq)
    mm,_,_ = nw.align(hs_seq,mm_seq)
    br,_,_ = nw.align(hs_seq,br_seq)
    tt,_,_ = nw.align(hs_seq,tt_seq)
    aligned = [gg,mm,br,tt]
    species = ['gallus', 'mus', 'balaeniceps', 'tursiops']

    print('Order of species from most to least similar:', [species[i] for i in np.argsort(aligned)][::-1])
    print('Alignment scores from most to least similar:', np.sort(aligned)[::-1])

if __name__ == "__main__":
    main()
