from Bio import SeqIO
DCT = {'A':'U','U':'A','C':'G','G':'C'}

def is_compl(base1: str, base2: str) -> bool:
    return DCT[base1] == base2

def cat(seq: str, memo: dict) -> int:
    l = len(seq)
    if (l < 2):
        return 1
    if seq in memo:
        return memo[seq]
    counter = 0
    for n in range(1, l, 2):
        if is_compl(seq[0], seq[n]):
            counter += cat(seq[1:n], memo) * cat(seq[n+1:], memo)
    memo[seq] = counter % 10**6
    return memo[seq]

def main():
    seq = SeqIO.read("rosalind_cat.txt", "fasta").seq
    memo = {}
    with open("out.txt", "w") as o:
        print(cat(seq, memo), file=o)

if __name__ == "__main__":
    main()