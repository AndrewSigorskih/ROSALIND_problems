import numpy as np
from Bio import SeqIO
from itertools import product

def EditDistanceAlignment(s1: str, s2:str) -> tuple:
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
        tab[i][j] = min(tab[i-1, j] +1,
                        tab[i, j-1] + 1,
                        tab[i-1, j-1] + cost)
    edit_distance = tab[m][n]
    # backtrack
    res1, res2 = "", ""
    i, j = m, n
    while ((i > 0)  and (j > 0)):
        left = tab[i-1, j]
        top = tab[i, j-1]
        topleft = tab[i-1, j-1]
        curmin = min(left, top, topleft)
        if (curmin == tab[i, j]): # no insertion or deletion
            res1 = s1[i-1] + res1
            res2 = s2[j-1] + res2
            i, j = i-1, j-1
        else:
            if (curmin == top) and ( curmin != left): # deletion
                res1 = "-" + res1
                res2 = s2[j-1] + res2
                j -= 1
            elif (curmin != top) and (curmin == left): # insertion
                res1 = s1[i-1] + res1
                res2 = '-' + res2
                i -= 1
            else: # substitution
                res1 = s1[i-1] + res1
                res2 = s2[j-1] + res2
                i, j = i-1, j-1 
    return (edit_distance, res1, res2)

def main():
    seq1, seq2 = (x.seq for x in SeqIO.parse("rosalind_edta.txt", "fasta"))
    res = EditDistanceAlignment(seq1, seq2)
    with open("out.txt", "w") as o:
        print(*res, sep='\n', file=o)

if __name__ == "__main__":
    main()