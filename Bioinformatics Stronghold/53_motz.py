from Bio import SeqIO
DCT = {'A':'U','U':'A','C':'G','G':'C'}

def is_compl(base1: str, base2: str) -> bool:
    return DCT[base1] == base2

def motz(seq: str, memo: dict) -> int:
    l = len(seq)
    if (l < 2):
        return 1
    if seq in memo:
        return memo[seq]
    counter = motz(seq[1:], memo)
    for k in range(1, l):  
        if is_compl(seq[0], seq[k]):
            counter += motz(seq[1:k], memo) * motz(seq[k+1:], memo)
    memo[seq] = counter % 10**6
    return memo[seq]
    
def main():
    seq = SeqIO.read("rosalind_motz.txt", "fasta").seq
    memo = {}
    with open("out.txt", "w") as o:
        print(motz(seq, memo), file=o)

if __name__ == "__main__":
    main()