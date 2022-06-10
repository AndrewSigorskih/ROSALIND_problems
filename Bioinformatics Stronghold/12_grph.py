from Bio import SeqIO
from itertools import product
Ovrlp = 3

records = list(SeqIO.parse("rosalind_grph.txt", "fasta"))
with open("out.txt", "w") as f:
    for i, j in product(range(len(records)), range(len(records))):
        if (i != j) and \
            (records[i].seq[-Ovrlp:] == records[j].seq[:Ovrlp]):
            print(records[i].id, records[j].id, file=f)
