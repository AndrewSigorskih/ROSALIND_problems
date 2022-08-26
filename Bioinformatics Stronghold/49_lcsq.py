from Bio import SeqIO
import numpy as np

import sys
# not needed if non-recursive baktrack is used
sys.setrecursionlimit(10_000)

# collect LCS
# https://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Reading_out_a_LCS
def Backtrack(L, X, Y, x, y):
    if (x == 0) or (y == 0):
        return ""
    if (X[x-1] == Y[y-1]):
        return Backtrack(L, X, Y, x-1, y-1) + X[x-1]
    if (L[x, y-1] > L[x-1, y]):
        return Backtrack(L, X, Y, x, y-1)
    return Backtrack(L, X, Y, x-1, y)

# non-recursive backtrack
def _Backtrack(L, X, Y, x, y):
    res = ""
    i, j = x, y
    while (i >0) and (j > 0):
        if X[i-1] == Y[j-1]:
            res += X[i-1]
            i -= 1
            j -= 1
        elif L[i-1, j] > L[i, j-1]:
            i -= 1
        else:
            j -= 1
    return res[::-1]

# Dynamic programming implementation of LCS problem
# Time Complexity: O(m*n)
def LongestCommonSubsequence(X, Y):
    m, n = len(X), len(Y)
    L = np.full((m+1, n+1), np.nan, dtype=np.int32)
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i, j] = 0
            elif X[i-1] == Y[j-1]:
                L[i, j] = L[i-1, j-1] + 1
            else:
                L[i, j] = max(L[i-1, j], L[i, j-1])
 
    return _Backtrack(L, X, Y, m, n)

def main():
    seq1, seq2 = (x.seq for x in SeqIO.parse("rosalind_lcsq.txt", "fasta"))
    with open("out.txt", "w") as o:
        print(LongestCommonSubsequence(seq1, seq2), file=o)

if __name__ == "__main__":
    main()