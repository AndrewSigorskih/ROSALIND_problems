# Dynamic programming - O(m*n) in time and space
# https://www.geeksforgeeks.org/shortest-possible-combination-two-strings/

import numpy as np
from itertools import product

def scsp(a: str, b: str) -> str:
    m = len(a)
    n = len(b)
    tab = np.zeros(shape=(m+1, n+1), dtype=int)
    # fill table
    for i,j in product(range(m+1), range(n+1)):
        if (not i):
            tab[i, j] = j
        elif (not j):
            tab[i, j] = i
        elif (a[i-1] == b[j-1]):
            tab[i, j] = 1 + tab[i-1, j-1]
        else:
            tab[i, j] = 1 + min(tab[i-1, j],
                                tab[i, j-1])
    # construct supersequence
    index = tab[m, n]
    res = ["" for _ in range(index)]
    i = m
    j = n
    # traverse table
    while (i > 0 and j > 0):
        if (a[i-1] == b[j-1]):
            res[index-1] = a[i-1]
            i -= 1
            j -= 1
            index -= 1
        elif (tab[i-1, j] < tab[i, j-1]):
            res[index-1] = a[i-1]
            i -= 1
            index -= 1
        else:
            res[index-1] = b[j-1]
            j -= 1
            index -= 1
    
    # Copy remaining characters of string a
    while (i > 0):
        res[index-1] = a[i-1]
        i -= 1
        index -= 1
 
    # Copy remaining characters of string b
    while (j > 0):
        res[index-1] = b[j-1]
        j -= 1
        index -= 1
 
    # return the result
    return "".join(res)

def main():
    with open("rosalind_scsp.txt", "r") as f:
        s1 = f.readline().strip('\n')
        s2 = f.readline().strip('\n')
    res = scsp(s1, s2)
    with open("out.txt", "w") as o:
        print(res, file=o)

if __name__ == "__main__":
    main()