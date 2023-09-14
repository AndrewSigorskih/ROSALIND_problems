from collections import Counter

with open ('rosalind_dna.txt', 'r') as infile:
    sequence = infile.readline().strip()
cnt = Counter(sequence)
print(cnt["A"], cnt["C"], cnt["G"], cnt["T"])
