import numpy as np
from Bio import SeqIO

with open("rosalind_bphr.txt", "r") as handle:
    q = int(handle.readline().strip())
    records = list(SeqIO.parse(handle, format="fastq"))
arr = np.array([x.letter_annotations["phred_quality"] for x in records])
print((arr.mean(axis=0) < q).sum())