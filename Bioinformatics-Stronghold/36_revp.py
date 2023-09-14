from Bio import SeqIO

def main():
    record = SeqIO.read("rosalind_revp.txt", "fasta").seq
    with open("out.txt", "w") as o:
        for k in range(4, 13):
            for i in range(len(record) - k + 1):
                if (str(record[i:i+k]) ==\
                     str(record[i:i+k].reverse_complement())):
                     print(i+1, k, file=o)

if __name__ == "__main__":
    main()