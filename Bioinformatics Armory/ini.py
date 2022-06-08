from Bio.Seq import Seq
with open ('rosalind_ini.txt', 'r') as infile:
    sequence = Seq(infile.readline().strip("\n"))
print(*[sequence.count(x) for x in "ACGT"], \
        sep=" ", end='\n')
