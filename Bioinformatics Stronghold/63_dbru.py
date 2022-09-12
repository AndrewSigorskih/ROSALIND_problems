from Bio.Seq import Seq

def revcmp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def get_edge(s: str) -> str:
    k = len(s)
    return f"({s[:k-1]}, {s[1:k]})"

def main():
    reads = set()
    with open("rosalind_dbru.txt", "r") as f:
        for line in f:
            seq = line.strip()
            reads.add(seq)
            reads.add(revcmp(seq))

    with open("out.txt", "w") as o:
        print(*sorted(get_edge(x) for x in reads), sep='\n', file=o)

if __name__ == "__main__":
    main()