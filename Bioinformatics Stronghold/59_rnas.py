DCT = {'A':('U'),
       'U':('A','G'),
       'C':('G'),
       'G':('C','U')
       }

def is_compl(base1: str, base2: str) -> bool:
    return base2 in DCT[base1]

def wobble_bonding(seq: str, memo: dict) -> int:
    length = len(seq)
    if (length <= 4):
        return 1
    if (seq in memo):
        return memo[seq]
    memo[seq] = wobble_bonding(seq[1:], memo)
    for pos in range(4, len(seq)):
        if is_compl(seq[0], seq[pos]):
            memo[seq] += wobble_bonding(seq[1:pos], memo) * wobble_bonding(seq[pos+1:], memo)
    return memo[seq]

def main():
    with open("rosalind_rnas.txt", "r") as f:
        seq = "".join((line.strip("\n") for line in f))
    memo = {}
    with open("out.txt", "w") as o:
        print(wobble_bonding(seq, memo), file=o)

if __name__ == "__main__":
    main()