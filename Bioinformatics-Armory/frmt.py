from Bio import Entrez, SeqIO

with open("rosalind_frmt.txt", "r") as f:
    ids = f.readline().strip("\n").replace(" ", ", ")

Entrez.email = "your_name@your_mail_server.com"
handle = Entrez.efetch(db="nucleotide", id=[ids], rettype="fasta")
records = list(SeqIO.parse(handle, "fasta")) 
with open("out_frmt.fasta", "w") as output_handle:
    SeqIO.write(min(records, key=len), output_handle, "fasta")