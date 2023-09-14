from Bio import SeqIO
def prefix(s):
    p = [0 for _ in range(len(s))]
    for i in range(1, len(s)):
        k = p[i - 1]
        while k > 0 and s[k] != s[i]:
            k = p[k - 1]
        if s[k] == s[i]:
            k += 1
        p[i] = k
    return p

def main():
    with open("rosalind_kmp.txt", "r") as f, \
        open("out.txt", "w") as o:
        s = SeqIO.read(f, "fasta").seq
        print(*prefix(s), file=o)

if __name__ == "__main__":
    main()