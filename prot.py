from Bio.Seq import Seq
with open('rosalind_prot.txt') as infile:
    matrix = Seq(infile.readline().strip('\n'))
with open('out.txt', 'w') as outfile:
    outfile.write(str(matrix.translate()).strip('*'))

