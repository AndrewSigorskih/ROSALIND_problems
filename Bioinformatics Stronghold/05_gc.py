from Bio import SeqIO
from Bio.SeqUtils import GC

records = SeqIO.parse("rosalind_gc.txt", format="fasta")
d = {record.name : GC(record.seq) for record in records}
k, v = max(d.items(), key=lambda x: x[1])
print(f"{k}\n{v:.6f}")