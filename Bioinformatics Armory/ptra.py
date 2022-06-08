from Bio import SeqIO, Seq

codes = Seq.CodonTable.ambiguous_generic_by_id.keys()
with open("rosalind_ptra.txt", "r") as f:
    dna, prot = (i.strip("\n") for i in f)

for N in codes:
    product = Seq.translate(dna, table=N, stop_symbol='')
    if prot in product:
        print(N)
        break