import numpy as np
from Bio import SeqIO

def hamming_dist(s1, s2):
    if (len(s1) != len(s2)):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(s1, s2)) / len(s1)

records = list(SeqIO.parse("rosalind_pdst.txt", format="fasta"))    
n = len(records)
arr = np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        arr[i,j] = hamming_dist(records[i].seq, records[j].seq)
        arr[j,i] = arr[i,j]

np.savetxt("out.txt", arr, fmt="%.5f") 