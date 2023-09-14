from math import factorial
from Bio import SeqIO

def main():
    record = SeqIO.read("rosalind_pmch.txt", "fasta")
    print(factorial(record.seq.count("A")) * \
        factorial(record.seq.count("C")))
        
if __name__ == "__main__":
    main()