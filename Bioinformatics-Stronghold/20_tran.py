from Bio import SeqIO

transitions, transversions = 0,0
s1, s2 = SeqIO.parse("rosalind_tran.txt", "fasta")
for c1, c2 in zip(s1.seq, s2.seq):
    if c1 != c2:
        if ({c1,c2} == {"A", "G"}) or ({c1, c2} == {"C", "T"}):
            transitions += 1
        else:
            transversions += 1
print(f"{(transitions / transversions):.11f}")