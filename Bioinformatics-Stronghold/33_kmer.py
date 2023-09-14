from itertools import product
from Bio import SeqIO
import re

alph = "ATGC"
K = 4

def main():
    record = SeqIO.read("rosalind_kmer.txt", "fasta")
    s = str(record.seq)
    kmers = sorted(["".join(x) for x in product(alph, repeat=K)])
    with open("out.txt", "w") as o:
        print(*[\
            sum(1 for _ in re.finditer(f'(?={pat})', s)) for pat in kmers\
                ], file=o)
        
if __name__ == "__main__":
    main()