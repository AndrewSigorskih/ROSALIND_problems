from Bio import SeqIO
mean = lambda x: sum(x)/len(x)

with open("rosalind_phre.txt", "r") as handle:
    thre = int(handle.readline())
    records = list(SeqIO.parse(handle, format="fastq"))
print(sum((1 for x in records \
    if mean(x.letter_annotations["phred_quality"]) < thre)))