from Newick import Tree

def main():
    with open("rosalind_sptd.txt", "r") as f:
        leaves, s1, s2 = [line.strip() for line in f]
    leaves = {leaf:i for i, leaf in enumerate(leaves.split())}
    t1, t2 = Tree(s1), Tree(s2)
    with open('out.txt', 'w') as o:
        print(Tree.rf_distance(t1, t2, leaves, unrooted=True), file=o)

if __name__ == "__main__":
    main()
