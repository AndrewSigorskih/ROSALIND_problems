# dynamic programming approach
import numpy as np
from Bio import SeqIO
from itertools import product

def edit_distance(s1: str, s2: str) -> int:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill table
    for i in range(m+1):
        tab[i, 0] = i
    for j in range(n+1):
        tab[0, j] = j
    for i, j in product(range(1, m+1), range(1, n+1)):
        if s1[i-1] == s2[j-1]:
            cost = 0
        else:
            cost = 1
        tab[i][j] = min(tab[i-1, j] + 1,
                        tab[i, j-1] + 1,
                        tab[i-1, j-1] + cost)
    return tab[m][n]
    
def main():
    seq1, seq2 = (x.seq for x in SeqIO.parse("rosalind_edit.txt", "fasta"))
    with open("out.txt", "w") as o:
        print(edit_distance(seq1, seq2), file=o)

if __name__ == "__main__":
    main()
