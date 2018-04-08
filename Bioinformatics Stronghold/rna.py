with open ('rosalind_rna.txt') as infile:
    sequence = infile.read()
print(sequence.replace("T", "U").strip())
