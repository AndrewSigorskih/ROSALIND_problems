import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from itertools import product

MAT = substitution_matrices.load("BLOSUM62")
MIN_INT = np.iinfo(np.int32).min
GAPOPEN = 11
GAPEXTEND = 1

def insert_gap(word, i):
    return word[:i] + "-" + word[i:]

def gaff(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    # here we will need to keep track 
    # of pointers for all 3 matrices
    # so 2 3*m*n tensors in total
    # lower: 0, middle:1, upper:2
    tab = np.zeros(shape=(3, m+1, n+1), dtype=int)
    pointers = np.zeros(shape=(3, m+1, n+1), dtype=int)
    # prepare table
    for i in range(1, m+1):
        tab[2, i, 0] = MIN_INT
        tab[1, i, 0] = -GAPOPEN
        tab[0, i, 0] = -GAPOPEN
    for j in range(1, n+1):
        tab[2, 0, j] = -GAPOPEN
        tab[1, 0, j] = -GAPOPEN
        tab[0, 0, j] = MIN_INT
    # fill table and keep backtrack
    for i, j in product(range(1, m+1), range(1, n+1)):
        # upper
        tmp = [tab[2, i, j-1] - GAPEXTEND, tab[1, i, j-1] - GAPOPEN]
        tab[2, i, j] = max(tmp)
        pointers[2, i, j] = tmp.index(tab[2, i, j])
        # lower
        tmp = [tab[0, i-1, j] - GAPEXTEND, tab[1, i-1, j] - GAPOPEN]
        tab[0, i, j] = max(tmp)
        pointers[0, i, j] = tmp.index(tab[0, i, j])
        # middle
        tmp = [tab[0, i, j], tab[1, i-1, j-1] + MAT[s1[i-1], s2[j-1]], tab[2, i, j]]
        tab[1, i, j] = max(tmp)
        pointers[1, i, j] = tmp.index(tab[1, i, j])

    # prepare indices and sequences for alignment
    i, j = m, n
    s1_al, s2_al = s1, s2
    scores = [tab[0, i, j], tab[1, i, j], tab[2, i, j]]
    max_score = max(scores)
    backtrack_pntr = scores.index(max_score)
    # backtrack
    while (i != 0) and (j != 0):
        if (backtrack_pntr == 0):
            # lower matrix
            if pointers[0, i, j] == 1:
                backtrack_pntr = 1
            i -= 1
            s2_al = insert_gap(s2_al, j)
        elif (backtrack_pntr == 1):
            # middle matrix
            if pointers[1, i, j] == 0:
                backtrack_pntr = 0
            elif pointers[1, i, j] == 2:
                backtrack_pntr = 2
            else:
                i -= 1
                j -= 1
        else:
            # upper matrix
            if pointers[2, i, j] == 1:
                backtrack_pntr = 1
            j -= 1
            s1_al = insert_gap(s1_al, i)
    # insert gaps in front if needed
    for _ in range(i):
        s1_al = insert_gap(s1_al, 0)
    for _ in range(j):
        s2_al = insert_gap(s2_al, 0)
    # return result
    return max_score, s1_al, s2_al


def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_gaff.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*gaff(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()