import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices
from itertools import product

MAT = substitution_matrices.load("BLOSUM62")
MIN_INT = np.iinfo(np.int32).min
GAPOPEN = 5

def gcon(s1: str, s2:str) -> int:
    m = len(s1)
    n = len(s2)
    # here we will need 3 matrices
    # just as with affine gap penalty, 
    # but gap extention is set to 0
    upper = np.zeros(shape=(m+1, n+1), dtype=int)
    middle = np.zeros(shape=(m+1, n+1), dtype=int)
    lower = np.zeros(shape=(m+1, n+1), dtype=int)
    for i in range(1, m+1):
        upper[i, 0] = MIN_INT
        middle[i, 0] = -GAPOPEN
        lower[i, 0] = -GAPOPEN
    for j in range(1, n+1):
        upper[0, j] = -GAPOPEN
        middle[0, j] = -GAPOPEN
        lower[0, j] = MIN_INT
    for i, j in product(range(1, m+1), range(1, n+1)):
        upper[i, j] = max(
            upper[i, j-1],
            middle[i, j-1] - GAPOPEN
        )
        lower[i, j] = max(
            lower[i-1, j],
            middle[i-1, j] - GAPOPEN
        )
        middle[i, j] = max(
            upper[i, j],
            lower[i, j],
            middle[i-1, j-1] + MAT[s1[i-1], s2[j-1]]
        )
    return middle[m, n]

def main():
    seq1, seq2 = (item.seq for item in SeqIO.parse("rosalind_gcon.txt", "fasta"))
    
    with open("out.txt", "w") as o:
        print(gcon(seq1, seq2), file=o)    

if __name__ == "__main__":
    main()