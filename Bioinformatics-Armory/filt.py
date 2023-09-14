from Bio import SeqIO

with open("rosalind_filt.txt", "r") as handle:
    q,p = map(int, handle.readline().split())
    records = list(SeqIO.parse(handle, format="fastq"))
p /= 100
count = 0
for record in records:
    if len(
        [x for x in record.letter_annotations["phred_quality"] if x >= q]
            ) / len(record) >= p:
        count += 1
print(count)