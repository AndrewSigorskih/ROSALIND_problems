import numpy as np
from Bio import SeqIO
from itertools import product

MOD = 2**27 -1

def ctea(s1: str, s2:str) -> int:
    m = len(s1)
    n = len(s2)
    # second table to keep track every step
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    counts = np.zeros(shape=(m+1, n+1), dtype=int)
    # init
    for i in range(m+1): # counts must start from 1 at [0, 0]
        tab[i, 0] = i
        counts[i, 0] = 1
    for j in range(1, n+1):
        tab[0, j] = j
        counts[0, j] = 1
    # fill tables
    for i, j in product(range(1, m+1), range(1, n+1)):
        cost = 1 if (s1[i-1] != s2[j-1]) else 0
        scores = [
            tab[i-1, j] + 1,
            tab[i, j-1] + 1,
            tab[i-1, j-1] + cost
        ]
        tab[i, j] = min(scores)
        # calc branches
        counts[i, j] += counts[i-1, j] if (tab[i, j] == scores[0]) else 0
        counts[i, j] += counts[i, j-1] if (tab[i, j] == scores[1]) else 0
        counts[i, j] += counts[i-1, j-1] if (tab[i, j] == scores[2]) else 0
        counts[i, j] = counts[i, j] % MOD
    return counts[m, n]

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_ctea.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(ctea(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()