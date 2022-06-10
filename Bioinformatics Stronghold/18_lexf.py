from itertools import product

with open("rosalind_lexf.txt", "r") as f:
    alph = f.readline().strip().split()
    k = int(f.readline().strip())
alph.sort()
with open("out.txt", "w") as f:
    for item in product(alph, repeat=k):
        print(*item, sep="", file=f)
