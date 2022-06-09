from Bio.Seq import Seq
with open ('rosalind_revc.txt', 'r') as infile:
    sequence = Seq(infile.readline().strip('\n'))
print(str(sequence.reverse_complement()))
