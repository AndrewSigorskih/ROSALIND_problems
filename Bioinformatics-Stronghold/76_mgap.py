import numpy as np
from Bio import SeqIO
from itertools import product

# this is basically the longest common subsequence problem:
# |LCS| of 2 sequences is the maximum number of aligned symbols
# the rest will be gaps in both sequences
# so, answer is |seq1| + |seq2| - 2*|LCS(seq1, seq2)|

def mgap(s1: str, s2:str) -> int:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    for i, j in product(range(1, m+1), range(1, n+1)):
        if (s1[i-1] == s2[j-1]):
            tab[i, j] = tab[i-1, j-1] +1
        else:
            tab[i, j] = max([
                tab[i-1, j],
                tab[i, j-1]
            ])

    return n + m - 2*tab[m, n]

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_mgap.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(mgap(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()