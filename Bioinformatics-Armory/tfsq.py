from Bio import SeqIO

records = SeqIO.parse("rosalind_tfsq.txt", "fastq")
SeqIO.write(records, "out_tfsq.fasta", "fasta")