from nis import match
import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from itertools import product

MAT = substitution_matrices.load("PAM250")
GAP_PEN = 5

def loca(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    pointers = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill tables
    for i, j in product(range(1, m+1), range(1, n+1)):
        scores = [
            tab[i-1, j] - GAP_PEN,
            tab[i, j-1] - GAP_PEN,
            tab[i-1, j-1] + MAT[s1[i-1], s2[j-1]],
            0 # end of local alignment
        ]
        tab[i, j] = max(scores)
        pointers[i, j] = scores.index(tab[i, j])
    # jump to highest scoring local alignment
    i, j = np.unravel_index(tab.argmax(), tab.shape)
    max_score = tab[i, j]
    s1_al, s2_al = s1[:i], s2[:j]
    # backtrack
    while (i != 0) and (j != 0): 
        if (pointers[i, j]) == 0:
            i -= 1
        elif (pointers[i, j]) == 1:
            j -= 1
        elif (pointers[i, j]) == 2:
            i -= 1
            j -= 1
        else: # (pointers[i, j] == 3) end of local aln
            break

    s1_al = s1_al[i:]
    s2_al = s2_al[j:]
    return max_score, s1_al, s2_al

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_loca.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*loca(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()