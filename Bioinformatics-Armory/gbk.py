from Bio import Entrez
with open("rosalind_gbk.txt", "r") as f:
    genus, start, end = (line.strip("\n") for line in f)
searchTerm = f'({genus}[Organism]) AND("{start}"[Publication Date]: "{end}"[Publication Date])'
Entrez.email = "sample@text.com"
handle = Entrez.esearch(db="nucleotide", term=searchTerm)
record = Entrez.read(handle)
print(record["Count"])