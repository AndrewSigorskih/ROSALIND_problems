from Bio.Seq import Seq

def revcmp(s: str) -> str:
    return str(Seq(s).reverse_complement())

def get_edge(s: str) -> tuple:
    k = len(s)
    return (s[:k-1], s[1:k])

# search for k that gives exactly two cycles
def gasm(reads: list) -> str:
    # try for all kmers of size len(read)...2
    for k_val in range(1, len(reads[0])):
        kmers = set()
        for read in reads:
            for i in range(k_val):
                kmers.add(read[i:len(read)+i-k_val+1])
                kmers.add(revcmp(read[i:len(read)-k_val+i+1]))
        edges = [get_edge(kmer) for kmer in kmers]
        assembled = []
        # try to assemble exactly 2 superstrings
        for _ in range(2): 
            cur_edge = edges.pop(0)
            result = cur_edge[0][-1]
            while cur_edge[1] in [item[0] for item in edges]:
                result += cur_edge[1][-1]
                [index] = [i for i, pair in enumerate(edges) if pair[0] == cur_edge[1]]
                cur_edge = edges.pop(index)
            assembled.append(result)
        if not edges: # empty -> all kmers are used -> exactly 2 cycles found
            return assembled[0]

def main():
    with open("rosalind_gasm.txt", "r") as f:
        reads = [line.strip() for line in f]
    res = gasm(reads)
    with open("out.txt", "w") as o:
        print(res, file=o)
if __name__ == "__main__":
    main()