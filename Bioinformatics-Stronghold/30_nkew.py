from Newick import TreeGraph

with open("rosalind_nkew.txt", "r") as f, open("out.txt", "w") as o:
    nwck, names, _ = (f.readline() for _ in range(3))
    while (nwck and names):
        nwck = nwck.strip()
        name1, name2 = names.strip().split()
        tree = TreeGraph(nwck)
        print(int(tree.get_distance(name1, name2)), end=" ", file=o)
        nwck, names, space = (f.readline() for _ in range(3))
    print("", end="\n", file=o)