import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from itertools import product

MAT = substitution_matrices.load("BLOSUM62")
MIN_INT = np.iinfo(np.int32).min
GAPOPEN = 11
GAPEXTEND = 1
# test datasets are too big, may not compute in time
# on some of them
def laff(s1: str, s2:str) -> tuple:
    m = len(s1)
    n = len(s2)
    max_score, max_i, max_j = MIN_INT, 0, 0
    upper = np.zeros(shape=(m+1, n+1), dtype=int)
    middle = np.zeros(shape=(m+1, n+1), dtype=int)
    lower = np.zeros(shape=(m+1, n+1), dtype=int)
    pntrs = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill tables
    for i, j in product(range(1, m+1), range(1, n+1)):
        upper[i, j] = max(
            upper[i, j-1] - GAPEXTEND,
            middle[i, j-1] - GAPOPEN
        )
        lower[i, j] = max(
            lower[i-1, j] - GAPEXTEND,
            middle[i-1, j] - GAPOPEN
        )
        scores = [
            lower[i, j],
            middle[i-1, j-1] + MAT[s1[i-1], s2[j-1]],
            upper[i, j],
            0 # index will be 3 -> backtrack stop
        ]
        middle[i, j] = max(scores)
        pntrs[i, j] = scores.index(middle[i, j])
        if (middle[i, j] > max_score):
            max_score = middle[i, j]
            max_i, max_j = i, j
    # backtrack
    i, j = max_i, max_j
    s1_al, s2_al = s1[:i], s2[:j]
    while (i > 0) and (j > 0):
        if (pntrs[i, j] == 0):
            i -= 1
        elif (pntrs[i, j] == 1):
            i -= 1
            j -= 1
        elif (pntrs[i, j] == 2):
            j -= 1
        else: # pntrs[i, j] == 3
            break
    s1_al, s2_al = s1_al[i:], s2_al[j:]
    return max_score, s1_al, s2_al

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_laff.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(*laff(seq1, seq2), sep='\n', file=o)    

if __name__ == "__main__":
    main()