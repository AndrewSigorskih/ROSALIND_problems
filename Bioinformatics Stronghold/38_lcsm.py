import numpy as np
from Bio import SeqIO
from itertools import product

#https://en.wikipedia.org/wiki/Longest_common_substring_problem
def LCSubstr(S : str, T : str):
    L = np.zeros((len(S), len(T)), dtype=int)
    z = 0
    ret = set()
    # Dynamic programming: O(|S|*|T|)
    for i, j in product(range(len(S)), range(len(T))):
        if S[i] == T[j]:
            if (i == 0) or (j == 0):
                L[i, j] = 1
            else:
                L[i, j] = L[i - 1, j - 1] + 1
            if L[i, j] > z:
                z = L[i, j]
                ret = {S[i-z+1 : i+1]} 
            elif L[i, j] == z:
                ret |= {S[i-z+1 : i+1]}
        else:
            L[i, j] = 0
    return ret

def updateMatches(S : str, ret: set):
    new_ret = set()
    for item in ret:
        new_ret |=  LCSubstr(S, item)
    maxlen = len(max(new_ret, key=len))
    new_ret = {x for x in new_ret if len(x) == maxlen}
    return new_ret

def main():
    records = (str(x.seq) for x in \
        SeqIO.parse("rosalind_lcsm.txt", "fasta"))
    ret = LCSubstr(next(records), next(records))
    for item in records:
        ret = updateMatches(item, ret)
    with open("out.txt", "w") as o:
        print(list(ret)[0], file=o)

if __name__ == "__main__":
    main()