with open ('rosalind_dna.txt') as infile:
    sequence = infile.read()
print(sequence.count("A"), sequence.count("C"), \
	sequence.count("G"), sequence.count("T"))

