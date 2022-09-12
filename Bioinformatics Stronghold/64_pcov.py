def get_edge(s: str) -> tuple:
    k = len(s)
    return (s[:k-1], s[1:k])

def main():
    reads = set()
    with open("rosalind_pcov.txt", "r") as f:
        for line in f:
            seq = line.strip()
            reads.add(seq) # no need for revc here

    edges = sorted(get_edge(x) for x in reads)
    # setup
    cur_edge = edges.pop(0)
    result = cur_edge[0][-1]
    while edges: # take out untill empty
        result += cur_edge[1][-1]
        [index] = [i for i, edge in enumerate(edges) if edge[0] == cur_edge[1]]
        cur_edge = edges.pop(index)

    with open("out.txt", "w") as o:
        print(result, file=o)

if __name__ == "__main__":
    main()