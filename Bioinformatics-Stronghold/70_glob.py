import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from itertools import product

MAT = substitution_matrices.load("BLOSUM62")
GAP_PEN = 5

def glob(s1: str, s2:str) -> int:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill table
    for i in range(m+1):
        tab[i, 0] = i * (-GAP_PEN)
    for j in range(n+1):
        tab[0, j] = j * (-GAP_PEN)
    for i, j in product(range(1, m+1), range(1, n+1)):
        tab[i,j] = max(
            tab[i-1, j] - GAP_PEN,
            tab[i, j-1] - GAP_PEN,
            tab[i-1, j-1] + MAT[s1[i-1], s2[j-1]]
        )
    return tab[m, n]

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_glob.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(glob(seq1, seq2), file=o)    

if __name__ == "__main__":
    main()