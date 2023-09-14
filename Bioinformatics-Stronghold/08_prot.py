from Bio import Seq

with open("rosalind_prot.txt", "r") as f:
    seq = Seq.Seq(f.readline().strip())
print(seq.translate(stop_symbol=''))