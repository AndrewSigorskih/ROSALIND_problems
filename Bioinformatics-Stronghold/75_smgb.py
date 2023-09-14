import numpy as np
from Bio import SeqIO
from itertools import product

GAP_PEN = 1

def smgb(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    pntrs = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill tables
    for i in range(m+1):
        pntrs[i, 0] = 2
    for j in range(n+1):
        pntrs[0, j] = 1
   
    for i, j in product(range(1, m+1), range(1, n+1)):
        scores = [
            tab[i-1, j-1] + (1 if s1[i-1] == s2[j-1] else -1),
            tab[i, j-1] - (GAP_PEN if i < m else 0),
            tab[i-1, j] - (GAP_PEN if j < n else 0)
        ]
        tab[i, j] = max(scores)
        pntrs[i, j] = scores.index(tab[i, j])

    # backtrack
    i, j = m, n
    s1_al, s2_al = "", ""
    max_score = tab[i, j]
    while (i > 0 ) or (j > 0):
        if (pntrs[i, j] == 0):
            s1_al += s1[i - 1]
            s2_al += s2[j - 1]
            i, j = i-1, j-1
        elif (pntrs[i, j] == 1):
            s2_al += s2[j - 1]
            s1_al += "-"
            j-= 1
        else: # pntrs[i, j] == 2
            s1_al += s1[i - 1]
            s2_al += "-"
            i -= 1

    return max_score, s1_al[::-1], s2_al[::-1]

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_smgb.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*smgb(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()