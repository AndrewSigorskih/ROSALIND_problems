from Bio import SeqIO

def main():
    targ, pat = map(lambda x: str(x.seq), \
        SeqIO.parse("rosalind_sseq.txt", "fasta"))
    idx = 0
    with open("out.txt", "w") as o:
        for i, c in enumerate(pat):
            if (i > 0) and (c == pat[i-1]):
                idx += 1
            idx = targ.find(c, idx)
            print(idx + 1, end = " ", file=o)
        print("", file=o)

if __name__ == "__main__":
    main()