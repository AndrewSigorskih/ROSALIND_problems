import numpy as np
from Bio import SeqIO
from itertools import product

def sims(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    pntrs = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill tables
    for i in range(m+1):
        tab[i, 0] = 0
        pntrs[i, 0] = 2
    for j in range(n+1):
        tab[0, j] = -j
        pntrs[0, j] = 1
    for i, j in product(range(1, m+1), range(1, n+1)):
        scores = [
            tab[i-1, j-1] + (1 if s1[i-1] == s2[j-1] else -1),
            tab[i, j-1] -1,
            tab[i-1, j] -1
        ]
        tab[i, j] = max(scores)
        pntrs[i, j] = scores.index(tab[i, j])
    # backtrack
    last_row = tab[:, n]
    max_score = np.max(last_row)
    i = np.argmax(last_row)
    j = n
    s1_al, s2_al = "", ""
    while (i > 0 ) and (j > 0):
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
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_sims.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*sims(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()