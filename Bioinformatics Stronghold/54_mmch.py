from Bio import SeqIO
from math import factorial
from collections import Counter

prob = lambda n, k: factorial(n) // factorial(n-k)

def mmch(seq: str) -> int:
    cnt = Counter(seq)
    return prob(max(cnt["A"], cnt["U"]), min(cnt["A"], cnt["U"])) *\
         prob(max(cnt["G"], cnt["C"]), min(cnt["G"], cnt["C"]))

def main():
    seq = SeqIO.read("rosalind_mmch.txt", "fasta").seq
    with open("out.txt", "w") as o:
        print(mmch(seq), file=o)

if __name__ == "__main__":
    main()