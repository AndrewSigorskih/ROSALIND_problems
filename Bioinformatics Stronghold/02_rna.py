with open ('rosalind_rna.txt', 'r') as infile:
    sequence = infile.readline().strip()
print(sequence.replace("T", "U"))
