from Bio import SeqIO

print(len([1 for x in SeqIO.parse("rosalind_rvco.txt", "fasta") \
        if x.seq == x.seq.reverse_complement()]))