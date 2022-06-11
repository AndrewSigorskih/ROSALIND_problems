# unnecessarily complicated solution
# for those who wonder how to convert 
# Bio.Phylo objects to nx Graphs
from Bio import Phylo
from io import StringIO
import networkx
# prepare tree for networkx conversion
def clade_names_fix(tree):
    for idx, clade in enumerate(tree.find_clades()):
        if not clade.name:
            clade.name=str(idx)
# fix due to nx accepting any object as node label
def node_names_fix(graph):
    mapper = {x:x.name for x in \
        list(networkx.dfs_preorder_nodes(graph))}
    networkx.relabel_nodes(graph, mapper, copy=False)

def get_dist(newick, names):
    handle = StringIO(newick)
    tree = Phylo.read(handle, "newick")
    clade_names_fix(tree)
    G = Phylo.to_networkx(tree)
    node_names_fix(G)
    spl = dict(networkx.all_pairs_shortest_path_length(G))
    return spl[names[0]][names[1]]

with open("rosalind_nwck.txt", "r") as f, open("out.txt", "w") as o:
    nwck, names, space = (f.readline() for _ in range(3))
    while (nwck and names):
        nwck = nwck.strip()
        names = names.strip().split()
        print(get_dist(nwck, names), sep="", end=" ", file=o)
        nwck, names, space = (f.readline() for _ in range(3))
    print("", end="\n", file=o)