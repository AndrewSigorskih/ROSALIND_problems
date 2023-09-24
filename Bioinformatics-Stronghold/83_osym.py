import numpy as np
from Bio import SeqIO
from itertools import product

GAP_PEN = 1

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
            tab[i-1, j-1] + (1 if (s1[i-1]==s2[j-1]) else -1)
        )
    return tab

def osym(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    forw_mat = glob(s1, s2)
    rev_mat = glob(s1[::-1], s2[::-1])
    tab = np.zeros(shape=(m, n), dtype=int)
    for i, j in product(range(m), range(n)):
        tab[i, j] = forw_mat[i, j] +\
                    (1 if (s1[i]==s2[j]) else -1) +\
                    rev_mat[m-i-1, n-j-1]
    return tab.max(), tab.sum()

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_osym.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*osym(seq1, seq2), file=o, sep='\n')    

if __name__ == "__main__":
    main()