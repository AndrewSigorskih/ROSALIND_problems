from Bio import Phylo
from io import StringIO
# in case not all distances present
def clade_attr_fix(tree):
    for idx, clade in enumerate(tree.find_clades()):
        if not clade.name:
            clade.name=str(idx)
        if not clade.branch_length:
            clade.branch_length = 0.0
        
def get_dist(newick, names):
    handle = StringIO(newick)
    tree = Phylo.read(handle, "newick")
    clade_attr_fix(tree)
    return tree.distance(names[0], names[1])

with open("rosalind_nkew.txt", "r") as f, open("out.txt", "w") as o:
    nwck, names, space = (f.readline() for _ in range(3))
    while (nwck and names):
        nwck = nwck.strip()
        names = names.strip().split()
        print(int(get_dist(nwck, names)), sep="", end=" ", file=o)
        nwck, names, space = (f.readline() for _ in range(3))
    print("", end="\n", file=o)