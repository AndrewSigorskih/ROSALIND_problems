import difflib
from Bio import SeqIO
from itertools import product

# doesnt work properly on large cases for some reason
def _get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    return s1[pos_a:pos_a+size]

def get_overlap(s1, s2):
    i = len(s2)
    while i > 0:
        if s2[:i] == s1[-i:]:
            return s1[-i:]
        if s2[-i:] ==  s1[:i]:
            return s1[:i]
        i -= 1
    return ""

def get_edge_list(records):
    res = {}
    for i,j in product(range(len(records)), repeat=2):
        if i == j:
            continue
        s1 = str(records[i].seq)
        s2 = str(records[j].seq)
        ovrlp = get_overlap(s1, s2)
        if len(ovrlp) >= min(len(s1)/2 -2, len(s2)/2 -2):
            if s1.startswith(ovrlp):
                res[f"{records[j].name}:{records[i].name}"] = ovrlp
            else:
                res[f"{records[i].name}:{records[j].name}"] = ovrlp
    return res

def get_path(edges):
    starts = {x.split(":")[0]:x for x in edges.keys()}
    ends = [x.split(":")[1] for x in edges.keys()]
    skew = [x for x in starts.keys() if x not in ends]
    if len(skew) != 1:
        print(f"Warning: graph is not a straight path! {skew}")
    skew = skew[0]
    return skew, starts

def main():
    records = list(SeqIO.parse("rosalind_long.txt", "fasta"))
    edges = get_edge_list(records)
    records = {x.name : str(x.seq) for x in records}
    start, path = get_path(edges)
    result = records[start]
    while True:
        edge = path.get(start, False)
        if not edge:
            break
        next = edge.split(":")[1]
        bridge = edges[edge]
        result += records[next].replace(bridge, "")
        start = next
    with open("out.txt", "w") as o:
        print(result, file=o)

if __name__ == "__main__":
    main()